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

subroutine read_ml_file()

#ifdef MLD_NDM
   use gen_com_m, ONLY: fnam,  rangml,lenfnam
#else
   use ondm_gen_com_m, ONLY: fnam,  rangml
#endif


   use module_afs, only: n_rbf_afs, n_rbf, n_cheb, afs_type, afs_type_bartok
   use ml_in_ndm_module, ONLY: ml_type, one_pi, &
      descriptor_type, &
      nd_data, nd_fingerprint, &
      iread_ml, isave_ml, &
      toy_model, debug,  &
      n_g2_eta, n_g2_rs, n_g3_eta, n_g3_zeta, n_g3_lambda, &
      l_max, j_max,  &
      alpha_soap, n_soap, lsoap, lsoap_fcut_wes, lsoap_diag, &
      weighted_auto, &
      seed, &
      write_desc, write_desc_dump, read_desc_dump, &
      weighted, weighted_3ch, &
      marginal_likelihood, strict_behler, descriptor_behler, &
      descriptor_g2, descriptor_g3, descriptor_afs, descriptor_g2_afs, descriptor_behler, descriptor_pow_so3, &
      descriptor_pow_so4, descriptor_pow_so3_3body, descriptor_bispectrum_so4, descriptor_g2_bispectrum_so4, descriptor_mtp, &
      descriptor_body, descriptor_ftnbody, descriptor_ace, &
      descriptor_g2_pow_so4, descriptor_tbind, descriptor_body, descriptor_zetabody, & 
      sign_stress, sign_stress_big_box, mld_fit_type, online_fit_batch_size, snap_class_constraints, &
      train_only, &
      mtp_poly_min, mtp_poly_max,   &
      desc_forces,  &
      eta_max_g2, eta_min_g2, rs_min_g2, rs_max_g2, &
      lsoap_diag, lsoap_norm, lsoap_lnorm, nspecies_soap, r_cut_width_soap, &
      lambda_krr, min_lambda_krr, max_lambda_krr, n_values_lambda_krr, ml_type_descriptors, &
      chemical_elements_invisible, mld_order, &
      mld_linear, mld_quadratic, mld_linear_extended, mld_polyc, polyc_n_poly, polyc_n_hermite, ml_type_analysis, &
      n_pca, classes_for_mcd, classes_train_for_sigma, rmat_dim, descriptor_milady, &
      mld_regularization_type, mld_regularization_type_home, &
      mld_type_quadratic, mld_type_quadratic_qml, mld_type_quadratic_qnml, &
      mld_type_quadratic_bilinear,  mld_type_quadratic_zaxa, mld_type_quadratic_ZX, &
      power_line, power_coeff_renorm, &
      write_design_matrix, write_test_design_matrix, mld_kernel, ml_type_krr, ml_type_basis, ml_type_nlinear, svd_rcond, &
      type_of_loss, lmask, mask_file, desc_file_format, dim_fourier_nbody, &
      delta_fourier_nbody, length_fourier_nbody, &
      Nfix, delta_fix_N_rcut, discrete_fix_N_rcut, fix_Nmax_neigh, hdf_positions, &
      fit_als
   use module_input_als_fit, only: als_nsteps, als_tol, als_ridge_k, als_rho, als_alpha_method, als_precond_type, &
      als_nnls_alpha, als_nnls_mode, als_nu_max, als_block_partition
   use module_body_desc, only: bond_dist_transform, bond_beta, l_body_order, body_D_max, bond_dist_ann
   use module_kernel, only: kernel_type, kernel_se, kernel_mc, kernel_ou, kernel_po, kernel_po_scaled, kernel_maha, &
      length_kernel, sigma_kernel,  write_kernel_matrix, &
      kernel_dump, kernel_dump_by_mahalanobis_norm,  &
      np_kernel_ref, np_kernel_full, np_omega, kernel_power, power_mcd, kernel_random, kernel_random_maha, &
      kernel_random_po, krff_type

   use module_bispectrum_so4, only: inv_r0_input, lbso4_diag
   use module_nlinear, only: order_nlinear, lnlinear_precond
   use k_cross_validation, only: kcross, n_kcross
   use def_kernels, ONLY: length_kse, sigma_kse
   use temporary_data_cov, only: dim_data
   use module_lbfgs_input, only: lbfgs_xtol, lbfgs_eps, &
      lbfgs_m_hess, lbfgs_max_steps, &
      lbfgs_print, lbfgs_gtol
   use module_optimization, only: lambda_krr_fake
   use time_check_general, only: debug_time, it_counter, time_read_db
   use module_db_poscar, only: iread_energy, ref_energy_per_element, procs_per_file
   use module_so3, only: n_rbf_so3, radial_pow_so3
   use module_neigh_local, only: r_cut, r_cut_width,  r_cut_in, r_cut_width_in,  type_fcut
   use module_ml_scalapack, only: scalapack_driver, debug_scalapack, nbr_predefined, nbc_predefined
   use module_zbl, only : zbl_potential, zbl_type, zbl_mode_default_k2b, zbl_mode_alone, r1_zbl, r2_zbl, rr_k2b
   use module_pair_cutoffs, only: setup_pair_cutoffs
   use module_variance, only : energy_variance, force_variance, stress_variance
   use module_cur, only : cur_eps, cur_kval, cur_cval, cur_rval
   use module_Sigma_matrix, only : method_sigma, method_sigma_by_Phi
   use module_optimization, only: optimize_weights_db, optimize_weights_chem, optimize_weights_L1, optimize_weights_L2, optimize_weights_Le, &
      optimize_ga_population, class_no_optimize_weights, no_class_weights, max_iter_optimize_weights, &
      factor_force_error, factor_energy_error, factor_stress_error, itopt
   use module_kernel_2b, only : activate_k2b, sigma_2b, delta_2b, np_radial_2b, r_cut_2b, r_cut_width_2b
   use module_chemical_species, only: fix_no_of_elements, chemical_elements,weight_per_element, weight_per_element_3ch, &
      img_num_ch, img_weighted, ws1, ws2, ws3, ws4
   use mld_mpi
   use mld_logger
  use module_ftnbody, only : dim_rff, length_rff, delta_rff, find_best_length, &
                             r_cut_ft2b, r_cut_width_ft2b, r_cut_ft3b, r_cut_width_ft3b, &
                             r_cut_ft4b, r_cut_width_ft4b, r_cut_ft5b, r_cut_width_ft5b, &
                             ftnbody_chem_mode, ftnbody_chem_rank, ftnbody_hash_channels, &
                             ftnbody_chem_rank_n, ftnbody_hash_channels_n, &
                             ftnbody_model, ftnbody_model_id, &
                             spip_degree_nbody, spip_degree_n, &
                             FTNBODY_MODEL_POLY, FTNBODY_MODEL_GRAMM, &
                             FTNBODY_MODEL_CPIP, FTNBODY_MODEL_SPIP
   use data_type, only: icall
   use module_db_setup, only: selection_type, selection_type_first, selection_type_first_start, &
                              selection_type_last, selection_type_random, &
                              mdb_file => db_file, mdb_path => db_path, drop_short_dist
   use module_ace_desc, only: l_ace_order, NU_LIMIT_MAX, ace_numax, ace_nmax_list, ace_kmax_list, ace_nmax, ace_lmax, ace_kmax, ace_lmax_list, & 
                              ace_radial_poly, ace_lambda, ace_lambda_list, zetaace_order, &
                              dim_delta_zetaace, MAX_ZETAACE_ORDER, l_ace_set_rcut, &
                              ace_rcut_in_list, ace_rcut_out_list, ace_rcut_width_in_list, ace_rcut_width_out_list, ace_gencg, & 
                              ace_chem, ace_radial_chem, ace_npoints_spline, &
                              ace_chem_low_rank, ace_chem_low_rank_q, ace_chem_low_rank_niter, ace_chem_low_rank_lambda, &
                              ace_svd_randomized, ace_svd_randomized_oversample, ace_svd_randomized_power_iter

   use module_kernel_zetabody, only: dim_grid_zetabody, dim_grid_zetabody, dim_length_zetabody, zetabody_order, & 
                                     r_cut_z2b, r_cut_z3b, r_cut_width_z2b, r_cut_width_z3b, dim_delta_zetabody, & 
                                     MAX_ZETABODY_ORDER                          
                             
   use module_tbind, only: tb_ham_ss, tb_ham_pp, tb_ham_sp, tb_ham_dd, &
                           tb_nn_max, tb_kmax_list, tb_nmax_list, tb_lmax_list, tb_lambda_list, &
                           tb_kmax_grid, tb_lmax_grid, &
                           tb_rcut_in_list, tb_rcut_out_list, tb_rcut_width_in_list, tb_rcut_width_out_list, &
                           tb_rcut_in_ij_list, tb_rcut_out_ij_list, &
                           tb_rcut_width_in_ij_list, tb_rcut_width_out_ij_list, &
                           tb_model_lambda, tb_model_trace, tb_model_equivb, tbind_model, &
                           tb_power_trace_list, &
                           tb_filter_type, tb_filter_order, tb_filter_ncenter, &
                           tb_filter_width, tb_filter_width_ratio, &
                           tb_filter_use_first_moment, tb_filter_use_real_part, tb_filter_use_imag_part, &
                           tb_filter_degree, tb_filter_spectral_padding, &
                           tb_filter_lmin_list, tb_filter_lmax_list, &
                           tb_nrf, tb_nrg, tb_g_type, tb_svd

   use mld_hdf5, only: db_hdf5
   use module_extxyz, only: db_xyz, xyz_detect_mode
   use module_json, only: db_json, json_detect_mode
   use module_pair_cutoffs, only: r_cut_in_pair_force, r1_zbl_pair_force, r2_zbl_pair_force, &
                                  setup_pair_cutoffs, write_pair_cutoffs_xml
   
   implicit none

   integer :: snap_type_quadratic, snap_type_quadratic_copy, mld_type_quadratic_copy, mld_type_quadratic_default
   integer :: snap_regularization_type, snap_regularization_type_copy, mld_regularization_type_copy, mld_regularization_type_default
   integer :: snap_fit_type, snap_fit_type_copy, mld_fit_type_copy, mld_fit_type_default
   integer :: snap_order, snap_order_copy, mld_order_copy, mld_order_default
   integer :: type_of_eml_file, type_of_eml_file_copy, desc_file_format_copy, desc_file_format_default
   logical :: optimize_weights_db_copy, optimize_weights_db_default, &
      optimize_weights, optimize_weights_copy, optimize_weights_default
   ! Local fixed-length db_path/db_file for namelist read (gfortran
   ! truncates character(:),allocatable in namelists). Copied to module
   ! variables (module_db_setup) after the read.
   character(len=512) :: db_path, db_file
   
   namelist /input_ml/ ml_type, &
      descriptor_type,  &
      nd_data, dim_data, nd_fingerprint, &
      kcross, n_kcross, &
      iread_ml, isave_ml, kernel_type, &
      marginal_likelihood, &
      toy_model, debug, &
      length_kse, sigma_kse, &
      n_g2_eta, n_g2_rs, n_g3_eta, n_g3_zeta, n_g3_lambda, &
      n_rbf, n_rbf_afs, n_rbf_so3, n_cheb, l_max, j_max, r_cut, &
      alpha_soap, n_soap, lsoap, lsoap_fcut_wes, &
      weighted_auto,  &
      seed, selection_type, &
      write_desc, write_desc_dump, read_desc_dump, &
      weighted, weighted_3ch, strict_behler, &
      optimize_weights, optimize_weights_db, optimize_weights_chem, class_no_optimize_weights, optimize_weights_L1, optimize_weights_L2, optimize_weights_Le, &
      optimize_ga_population, &
      factor_force_error, factor_energy_error, factor_stress_error, &
      sign_stress, sign_stress_big_box, &
      procs_per_file, &
   ! this is in double for the moment will be removed ...
      mld_fit_type, snap_fit_type, online_fit_batch_size, &
      snap_class_constraints, train_only, mtp_poly_min, mtp_poly_max, lbso4_diag, max_iter_optimize_weights, inv_r0_input, &
      desc_forces, &
      eta_max_g2, eta_min_g2, rs_min_g2, rs_max_g2, &
      nspecies_soap,  lsoap_diag, lsoap_norm, lsoap_lnorm, &
      lambda_krr, min_lambda_krr, max_lambda_krr, n_values_lambda_krr, &
      fix_no_of_elements, weight_per_element, weight_per_element_3ch, chemical_elements, chemical_elements_invisible, &
   ! this is in double for the moment will be removed ...
      snap_order, mld_order, &
      polyc_n_poly, polyc_n_hermite, n_pca, &
      classes_for_mcd, rmat_dim, afs_type, dim_fourier_nbody, length_fourier_nbody, delta_fourier_nbody, &
      ftnbody_chem_mode, ftnbody_chem_rank, ftnbody_hash_channels, ftnbody_model, spip_degree_nbody, &
   ! this is in double for the moment will be removed ...
      mld_regularization_type, snap_regularization_type, &
   ! this is in double for the moment will be removed ...
      mld_type_quadratic, snap_type_quadratic, &
      scalapack_driver,  &
      power_line, power_coeff_renorm, &
      body_D_max, l_body_order, bond_dist_transform, bond_beta, bond_dist_ann, &
      write_design_matrix,  write_test_design_matrix, write_kernel_matrix, &
      length_kernel, sigma_kernel, kernel_dump, np_kernel_ref, np_kernel_full, np_omega, &
      kernel_power, order_nlinear, debug_time, &
      lbfgs_xtol, lbfgs_eps, lbfgs_m_hess, lbfgs_max_steps, lbfgs_print, lbfgs_gtol, power_mcd, &
      svd_rcond, iread_energy, radial_pow_so3, ref_energy_per_element, &
      r_cut_width_soap, type_of_loss, &
      nbc_predefined, nbr_predefined, lmask, mask_file, &
   ! this is in double for the moment will be removed ...
      type_of_eml_file, desc_file_format, &
      zbl_potential, zbl_type, r1_zbl, r2_zbl, method_sigma, cur_eps, cur_kval, cur_cval, cur_rval, &
      debug_scalapack, energy_variance, force_variance, stress_variance, &
      activate_k2b, sigma_2b, delta_2b, np_radial_2b, type_fcut, r_cut_width, r_cut_in, r_cut_width_in, &
      r_cut_2b, r_cut_width_2b, img_num_ch, img_weighted, ws1, ws2, ws3, ws4, krff_type, &
      r_cut_ft2b, r_cut_width_ft2b, r_cut_ft3b, r_cut_width_ft3b, &
      r_cut_ft4b, r_cut_width_ft4b, r_cut_ft5b, r_cut_width_ft5b, find_best_length, &
      Nfix, delta_fix_N_rcut, discrete_fix_N_rcut, fix_Nmax_neigh, & 
      l_ace_order, ace_numax, ace_nmax_list, ace_lmax_list, ace_kmax_list, ace_lambda_list, ace_radial_poly, & 
      dim_grid_zetabody,   dim_grid_zetabody, dim_length_zetabody , zetabody_order, & 
      r_cut_z2b, r_cut_z3b, r_cut_width_z2b, r_cut_width_z3b, dim_delta_zetabody, &
      zetaace_order, dim_delta_zetaace, & 
      l_ace_set_rcut,  &
      ace_rcut_in_list, ace_rcut_out_list, ace_rcut_width_in_list, ace_rcut_width_out_list, &
      ace_gencg, ace_chem, ace_radial_chem, ace_npoints_spline, &
      ace_chem_low_rank, ace_chem_low_rank_q, ace_chem_low_rank_niter, ace_chem_low_rank_lambda, &
      ace_svd_randomized, ace_svd_randomized_oversample, ace_svd_randomized_power_iter, &
      tb_ham_ss, tb_ham_pp, tb_ham_sp, tb_ham_dd, &
      tb_nn_max, tb_kmax_list, tb_nmax_list, tb_lmax_list, tb_lambda_list, &
      tb_kmax_grid, tb_lmax_grid, &
      tb_rcut_in_list, tb_rcut_out_list, tb_rcut_width_in_list, tb_rcut_width_out_list, &
      tb_rcut_in_ij_list, tb_rcut_out_ij_list, tb_rcut_width_in_ij_list, tb_rcut_width_out_ij_list, &
      tb_model_lambda, tb_model_trace, tb_model_equivb, tbind_model, tb_power_trace_list, &
      tb_filter_type, tb_filter_order, tb_filter_ncenter, tb_filter_width, tb_filter_width_ratio, &
      tb_filter_use_first_moment, tb_filter_use_real_part, tb_filter_use_imag_part, &
      tb_filter_degree, tb_filter_spectral_padding, tb_filter_lmin_list, tb_filter_lmax_list, &
      tb_nrf, tb_nrg, tb_g_type, tb_svd, &
      db_path, db_file, db_hdf5, hdf_positions, drop_short_dist, &
      als_nsteps, als_tol, als_ridge_k, als_rho, als_alpha_method, als_precond_type, &
      als_nnls_alpha, als_nnls_mode, als_nu_max, &
      r_cut_in_pair_force, r1_zbl_pair_force, r2_zbl_pair_force


   character(len=80)  :: text, append_infos
   character(len=:), allocatable    :: tmp, fnamtin
   integer  :: luml, none_class, i, nitems, nitems2, int_local, ntest, ierr
   logical  :: lexist 

   time_read_db = 0.0d0 
   zetabody_order = 1 
   !
   dim_grid_zetabody(2) = " 30 "
   dim_grid_zetabody(3)=" 10 10 20 " 
   !
   dim_length_zetabody(2) = " 0.5d0 "
   dim_length_zetabody(3) = " 0.5d0 0.5d0 0.5d0 "
   !
   r_cut_z2b = -5.3d0       
   r_cut_width_z2b = -1.5d0
   r_cut_z3b = -5.3d0
   r_cut_width_z3b = -1.5d0
   !
   ! should have the dimension of zeta_order: 
   dim_delta_zetabody(2) = " 0.1d0 1.d0 10.d0"
   dim_delta_zetabody(3) = " 0.1d0 1.d0 10.d0"


   zetaace_order = 1 
   if (allocated(l_ace_order)) deallocate(l_ace_order) ; allocate(l_ace_order(NU_LIMIT_MAX))
   l_ace_order(1:NU_LIMIT_MAX) = .false.
   ace_numax = 3
   ace_nmax_list = " 5 3 1"
   ace_lmax_list = " 0 4 2"
   ace_kmax_list = " 3 3 3"
   ace_lambda_list = " 2.d0 2.d0 2.d0"
   ace_radial_poly = 2 ! 1 - powPftouny ; 2 - expPaftouny ; 3 - simpBessel 
   ace_chem = 1 ! other option 0 first version, 1- standard ACE version, 2 - TS version
   ace_radial_chem = 1 ! 1 - RALF poly version (default) 2 - POD TS VERSION  3 - HSVD  5 - HSVD random projection contraction
   ace_npoints_spline = 1000 ! number of points for the radial spline grid
   ace_chem_low_rank = 0       ! 0 = disabled, 1 = enable low-rank chemical compression
   ace_chem_low_rank_q = 4     ! compression rank Q
   ace_chem_low_rank_niter = 50 ! number of ALS iterations
   ace_chem_low_rank_lambda = 1.0d-6 ! ridge regularization for ALS
   ace_svd_randomized = 0             ! 0 = exact dgesvd (default), 1 = randomized truncated SVD for HSVD/BLOCK_HSVD radial (recommended when the number of species is large, e.g. > ~30-50)
   ace_svd_randomized_oversample = 10 ! oversampling added to kmax for the random sketch dimension
   ace_svd_randomized_power_iter = 2  ! power iterations to sharpen the singular spectrum of the sketch
   dim_delta_zetaace(:)  = " 0.1d0 1.d0 0.1d0"
   l_ace_set_rcut = .false. 
   ace_rcut_in_list = " 1.2d0 1.2d0 1.2d0 "  
   ace_rcut_out_list = " 5.d0 5.d0 5.d0 "  
   ace_rcut_width_in_list = " 0.4d0 0.4d0 0.4d0 "
   ace_rcut_width_out_list = " 0.5d0 0.5d0 0.5d0 "
   ace_gencg = 2 ! 1 - use the draft method , 2 - Dusson-Ortner SVD
   ! acekmax = 3

   ml_type = 0             ! -2 analysis -1 write descriptors  1 - KRR ; 0 SNAP_1 ; 2 - GP (not yet implemented).
   iread_ml = 0            ! 0 compute in the fly, 1 read from previous run
   isave_ml = 0            ! 0 do nothing, 1 write on the HDD and run, 2 MPI and threading ...
   db_file = "db_model.in"
   db_path = "DB/"
   drop_short_dist = -1.d0 ! <=0: keep all configs; >0: drop configs with an interatomic distance below this threshold (Angstrom)
#if(MLD_HDF5)   
   db_hdf5 = .false.       ! input files are or not in hdf5 format. 
#endif 
   nd_data = 0             ! the number of files in the repository
   dim_data = 100          ! the dimension of the database

   hdf_positions = .true. 
   procs_per_file = nb_procsml      ! the number of MPI process that will handle a single db file. Used to define MPI grid
   nbc_predefined = 32
   nbr_predefined = 32
   desc_forces = .true.
   strict_behler = .false.
   kcross = .false.
   marginal_likelihood = .false.
   n_kcross = 0
   toy_model = .false.
   classes_for_mcd = " 01 02"
   iread_energy = 2

   lambda_krr = -1.d0
   n_values_lambda_krr = 21
   min_lambda_krr = 1.d-10
   max_lambda_krr = 1.d+10

   !some kernel stuffs.
   kernel_type = 4
   length_kse = 1.d0       ! hasta will be eliminated for compatibility reason still there
   sigma_kse = 1.d0        ! hasta - to be eliminated
   kernel_dump = kernel_dump_by_mahalanobis_norm
   np_kernel_ref = 200
   np_kernel_full = 800
   np_omega = 4
   kernel_power = 2.d0
   length_kernel = sqrt(1.d0/2.d0)
   sigma_kernel = 0.d0
   power_mcd = 0.05d0     ! selection of the kernel
   krff_type = 1
   itopt=0
   icall = 0
   ! this is the default value ...
   svd_rcond = 100*epsilon(1.d0)
   svd_rcond = -1

   ! ALS-Ridge defaults (mld_fit_type = 5)
   als_nsteps = 10
   als_tol = -1.0d0
   als_ridge_k = 1.0d-3
   als_rho = 1.0d-6
   als_alpha_method = 1      ! 0 = fixed alpha, 1 = learning alpha
   als_precond_type = 0      ! 0 = Frobenius norm, 1 = SVD, 2 = flat (lambda=lambda_krr)
   als_nnls_alpha = .false.  ! enforce alpha >= 0 via NNLS
   als_nnls_mode = 0         ! 0 = Lawson-Hanson, 1 = BK warm-start
   als_nu_max = 1            ! number of blocks (must be set by user for fit_als)

#ifdef MLD_NDM
   fnamtin = fnam(1:lenfnam)//'.ml'
#else
   fnamtin = fnam//'.ml'
#endif
   descriptor_type = 1
   snap_class_constraints = "02"
   sign_stress = 1.d0
   sign_stress_big_box = 1.d0
   train_only = .false.
   optimize_ga_population = 40

   optimize_weights = .false.
   optimize_weights_default = .false.
   optimize_weights_db = .false.
   optimize_weights_db_default = .false.

   optimize_weights_chem = .false.
   optimize_weights_L2 = .false.
   optimize_weights_L1 = .false.
   optimize_weights_Le = .false.
   class_no_optimize_weights = " "
   max_iter_optimize_weights = 40
   factor_force_error = 1.d0
   factor_energy_error = 1.d0
   factor_stress_error = 1.d0
   n_pca = 3

   !g2
   n_g2_eta = 3
   n_g2_rs = 1

   !eta_max_g2 = 0.8d0
   rs_min_g2 = 0.d0
   rs_max_g2 = 0.d0
   eta_max_g2 = 200.d0
   eta_min_g2 = 0.2d0
   !g3
   n_g3_eta = 3
   n_g3_zeta = 2
   n_g3_lambda = 2

   !mtp
   mtp_poly_min = 2
   mtp_poly_max = 5


   !afs
   n_rbf = 4               ! afs, pow_so3, bispectrum_so3
   n_cheb = 5              ! afs
   afs_type = afs_type_bartok



   l_max = 4               ! pow_so3, bispectrum_so3, soap
   radial_pow_so3 = 1      ! 1   Bartok, 2  de Gironcolli

   j_max = 1.5             ! pow_so4, bispectrum_so4
   inv_r0_input = (1.d0 - 0.02d0/one_pi)
   lbso4_diag = .false.
   ! jj_max=int(2*j_max)
   r_cut = 5.d0            ! active g2, g3, rbf, pow_so4, bispectrum_so4, soap
   r_cut_width = 0.5d0
  r_cut_ft2b = 5.d0 
  r_cut_width_ft2b = 0.5d0
  r_cut_ft3b = 5.d0 
  r_cut_width_ft3b = 0.5d0
  r_cut_ft4b = 5.d0 
  r_cut_width_ft4b = 0.5d0
  r_cut_ft5b = 5.d0 
  r_cut_width_ft5b = 0.5d0

   r_cut_2b = -1.d0
   r_cut_width_2b = 0.5d0
   r_cut_width_soap = -0.5d0
   r_cut_in=-1.2d0
   r_cut_width_in= 0.4d0

   type_fcut = 2           ! 1, 2 or 3

   n_soap = 2              ! radial part of soap
   lsoap = .false.         ! soap
   lsoap_fcut_wes = .false.                         ! acd
   !debug purposes
   !alpha_soap=2d0      ! soap
   !atom_sigma=0.5, alpha_soap = 0.5/atom_sigma**2
   !atom_sigma_soap=0.5d0
   !alpha_soap = 0.5d0/(atom_sigma_soap**2)
   alpha_soap = 2.d0
   lsoap_diag = .false.    !
   lsoap_norm = .true.     !
   lsoap_lnorm = .true.
   nspecies_soap = 1       ! 1 - is our style, 2 - is Gabor style


   ! milady
   rmat_dim = 20
   power_line = 2.d0
   power_coeff_renorm = 0.333d0
   img_weighted = .false.
   img_num_ch = 1
   ws1= "1.d0"
   ws2= "1.d0"
   ws3= "1.d0"
   ws4= "1.d0"

   ! body developpement
   do i = 1, 5
      body_D_max(i) = i
   end do
   l_body_order(:) = .true.
   l_body_order(5) = .false.
   l_body_order(6) = .false.
   ! 1 - dist, 2 - exp 3 - inverse
   bond_dist_transform = 3
   bond_beta = 2.0
   bond_dist_ann = 1.0



   polyc_n_poly = 2
   polyc_n_hermite = 2
   selection_type = 3      ! 1 - selects first "ns" elements of the database
   ! 2 - selects last "ns" elements of the database
   ! 3 - selects randomly "ns" subsets of "kelem" elements of the database
   seed = 14792435               ! seed for random generator
   write_desc = .false.    ! stores ATOMIC descriptors
   write_desc_dump = .false.                        ! stores ATOMIC / OTHER INFORMATIONS for descriptors used for further utilisation
   read_desc_dump = .false.                         ! read  descriptors from ATOMIC / OTHER INFORMATIONS. No new calculation for descriptors are performed.
   write_design_matrix = .false.
   write_test_design_matrix = .false.

   write_kernel_matrix = .false.
   debug_time = .false.
   it_counter(:) = 0


   weighted = .false.      ! if true two channels version
   weighted_3ch = .false.  ! if true three channel version
   weighted_auto = .false. ! if true two channels version
   fix_no_of_elements = 1
   chemical_elements = "Fe"
   chemical_elements_invisible = " "
   weight_per_element = "1.d0"
   weight_per_element_3ch = "1.d0"
   ref_energy_per_element = "0.d0"
   scalapack_driver = .false.
   debug_scalapack = .false.
   nd_fingerprint = 3

   order_nlinear = 2
   lnlinear_precond = .true.


   lbfgs_xtol = 1.d-16
   lbfgs_eps = 1.d-2
   lbfgs_m_hess = 40
   lbfgs_max_steps = 20000
   lbfgs_print(1) = 100
   lbfgs_print(2) = 0
   lbfgs_gtol = 4.d-1

   !TODO ace_lambda = 2.d0
   type_of_loss = 1
   lmask = .false.
   mask_file=.false.

   zbl_potential = .false.
   zbl_type = zbl_mode_default_k2b
   r1_zbl = 1.d0
   r2_zbl= 2.2d0

   force_variance = .false.
   energy_variance = .false.
   stress_variance = .false.

   method_sigma = method_sigma_by_Phi

   cur_kval = -1
   cur_rval = -1
   cur_cval = -1
   cur_eps = 1.d0

   activate_k2b=.false.
   sigma_2b=0.05d0
   delta_2b=1.d0
   np_radial_2b=30

   !fixed number of neighbours
   Nfix=.false.
   discrete_fix_N_rcut=10
   fix_Nmax_neigh=8
   delta_fix_N_rcut=0.2d0

   ! things that I will change the name without people realise
   !mld_fit_type = 4
   !snap_class_constraints = "02"
   !snap_order
   !mld_kernel

   !snap_type_quadratic -> mld_type_quadratic
   mld_type_quadratic = -777
   snap_type_quadratic = -777
   mld_type_quadratic_copy  = mld_type_quadratic
   snap_type_quadratic_copy = snap_type_quadratic
   mld_type_quadratic_default = 1

   !snap_regularization_type -> mld_regularization_type
   mld_regularization_type = -777
   snap_regularization_type = -777
   mld_regularization_type_copy  = mld_regularization_type
   snap_regularization_type_copy = snap_regularization_type
   mld_regularization_type_default = 0

   ! sbap_fit_type -> mld_fit_type
   mld_fit_type = -777
   snap_fit_type = -777
   mld_fit_type_copy  = mld_fit_type
   snap_fit_type_copy = snap_fit_type
   mld_fit_type_default = 4

   online_fit_batch_size = 4 ! configs batched per collective scatter in the online Gram fit (mld_fit_type=6); 1 = original per-config behaviour

   mld_order = -777
   snap_order = -777
   mld_order_copy  = mld_order
   snap_order_copy = snap_order
   mld_order_default = 1

   desc_file_format= -777
   type_of_eml_file = -777
   desc_file_format_copy = desc_file_format
   type_of_eml_file_copy = type_of_eml_file
   desc_file_format_default = 1

   dim_fourier_nbody = " 0 0 2000 0 0"
   length_fourier_nbody = " 0.7d0 0.7d0 0.7d0 0.7d0 0.7d0 "
   delta_fourier_nbody = " 1.d0 1.d0 1.d0 1.d0 1.d0 "
   ftnbody_chem_mode = 0
   ftnbody_chem_rank = " 0 0 0 0 0 "
   ftnbody_hash_channels = " 0 0 0 0 0 "
   ftnbody_model = 'poly'
   spip_degree_nbody = " 0 4 4 4 3 "
   find_best_length = .true.


   inquire(file=fnamtin, exist=lexist)
   if (.not. lexist) then
      call log_critical('read_ml: input file not found: '//fnamtin)
      call log_critical('read_ml: please provide the input file')
      call mld_mpi_finalize('read_ml: please provide the input file')
      stop
   end if

   open (newunit=luml, file=fnamtin, status='old', iostat=ierr)
   if (ierr /= 0) then
      call log_critical('read_ml: error opening file: '//fnamtin)
      call log_critical('read_ml: check file permissions or existence')
      call mld_mpi_finalize('read_ml: check file permissions or existence')
      stop
   end if

   block
      character(len=512) :: errmsg
      errmsg = ' '
      read (luml, nml=input_ml, iostat=ierr, iomsg=errmsg)
      if (ierr /= 0) then
         call log_critical('read_ml: error reading namelist from file: '//fnamtin)
         if (ierr < 0) then
            call log_critical('read_ml: end of file reached (file might be empty)')
         else
            call log_critical('read_ml: format error in namelist input_ml')
         end if
         call log_critical('read_ml: iomsg = '//trim(errmsg))
         call mld_mpi_finalize('read_ml: fix the input file')
         stop
      end if
   end block
   close (luml)

   ! Copy local fixed-length vars to module allocatable variables
   ! (workaround for gfortran namelist truncation of character(:),allocatable)
   mdb_path = trim(db_path)
   mdb_file = trim(db_file)

   ! ---- Detect extended XYZ mode from db_path extension ----
   call xyz_detect_mode()
   if (db_xyz .and. rangml == 0) then
      write(6, '("ML: Extended XYZ mode detected from db_path = ", a)') trim(db_path)
   end if

   ! ---- Detect MPtrj-style JSON mode from db_path extension ----
   call json_detect_mode()
   if (db_json .and. rangml == 0) then
      write(6, '("ML: JSON dataset mode detected from db_path = ", a)') trim(db_path)
   end if

   !refix the double input ...

   !fix ... mld_type_quadratic
   call fix_double_input_for_integer(" mld_type_quadratic ", mld_type_quadratic, mld_type_quadratic_copy, &
      " snap_type_quadratic ", snap_type_quadratic, snap_type_quadratic_copy, &
      mld_type_quadratic_default)
   !fix ... mld_regularization_type
   call fix_double_input_for_integer(" mld_regularization_type ", mld_regularization_type, mld_regularization_type_copy, &
      " snap_regularization_type ", snap_regularization_type, snap_regularization_type_copy, &
      mld_regularization_type_default)
   !fix ... mld_fit_type
   call fix_double_input_for_integer(" mld_fit_type ", mld_fit_type, mld_fit_type_copy, &
      " snap_fit_type ", snap_fit_type, snap_fit_type_copy, &
      mld_fit_type_default)

   !fix ... mld_order
   call fix_double_input_for_integer(" mld_order ", mld_order, mld_order_copy, &
      " snap_order ", snap_order, snap_order_copy, &
      mld_order_default)

   !fix ... mld_order
   call fix_double_input_for_integer(" desc_file_format ", desc_file_format, desc_file_format_copy, &
      " type_of_eml_file ", type_of_eml_file, type_of_eml_file_copy, desc_file_format_default)

   !fix ... optimize_weights
   call fix_double_input_for_logical(" optimize_weights_db ", optimize_weights_db, optimize_weights_db_copy, &
      " optimize_weights ", optimize_weights, optimize_weights_copy, &
      optimize_weights_db_default)

   !end fix ..................
   ! TODO FIX THAT
   n_rbf_afs = n_rbf
   n_rbf_so3 = n_rbf
   ! TODO FIX THAT
   if (r_cut_width_soap < 0.d0 ) then
      r_cut_width_soap = r_cut_width
   end if

   if (mld_regularization_type == mld_regularization_type_home) call set_grid_lambda_krr

   ! ALS-Ridge: allocate block partition vector after reading als_nu_max
   ! Note: als_block_partition is NOT in the namelist (it is allocatable).
   ! The actual partition boundaries are resolved later in the dispatch code
   ! (parameters.F90) where the design matrix dimensions are known.
   if (mld_fit_type == fit_als) then
     if (als_nu_max < 1) als_nu_max = 1
     if (allocated(als_block_partition)) deallocate(als_block_partition)
     allocate(als_block_partition(2*als_nu_max))
     als_block_partition(:) = 0
   end if


   !if ((mld_fit_type==fit_home_hb) .and.(optimize_weights)) then
   if (zbl_potential) then
      if (zbl_type /= zbl_mode_default_k2b .and. zbl_type /= zbl_mode_alone) then
         call mld_mpi_abort('MLD fatal: zbl_type must be zbl_mode_default_k2b or zbl_mode_alone. Got: '//vtoa(zbl_type))
      end if

      if (r1_zbl == r2_zbl) then
         call mld_mpi_abort('MLD fatal: zbl_potential is true and the two r1_zbl and r2_zbl are equal. Should not.  '//vtoa(r1_zbl)//vtoa(r2_zbl))
      end if
      if (r1_zbl > r2_zbl) then
         call mld_mpi_abort('MLD fatal: zbl_potential is true and r1_zbl should be lower than r2_zbl. Is not the case. '//vtoa(r1_zbl)//vtoa(r2_zbl))
      end if

      if (r_cut_in <= 0 ) then
         ! Negative r_cut_in signals automatic pair-specific mode.
         ! setup_pair_cutoffs() will compute per-pair values from covalent radii.
         call log_info('read_ml: r_cut_in is negative ('//vtoa(r_cut_in)//') => automatic pair-specific mode'// &
            ' (will be computed from covalent radii in setup_pair_cutoffs)')
      end if

      if (r_cut_in > 0 ) then
         if (r_cut_in <= r1_zbl) then
           call log_critical('MLD fatal: zbl_potential is true and the r1_zbl is larger than r_cut_in. Should not. ')
           call log_critical('For the moment r1_zbl and r_cut_in have the values '//vtoa(r1_zbl)//vtoa(r_cut_in)//'  ') 
           call mld_mpi_abort('Change the value of r_cut_in such that to be larger than r1_zbl and lower than r2_zbl')
         end if
         if (r_cut_in >= r2_zbl) then
            call log_critical('MLD fatal: zbl_potential is true and the r_cut_in is larger than r2_zbl. Should not. ')
            call log_critical('For the moment r2_zbl and r_cut_in have the values '//vtoa(r2_zbl)//vtoa(r_cut_in)//'  ') 
            call mld_mpi_abort('Change the value of r_cut_in such a way to be larger than r1_zbl and lower than r2_zbl')
         end if
      end if

      if (zbl_type == zbl_mode_default_k2b) then
         ! zbl_mode_default_k2b: k2b bridge mode (original behavior)
         rr_k2b = r_cut_in

         if (.not.activate_k2b) then
            call log_warning('The initial kernel 2b was descativated and will be activated. Why ? Why not ! ')

            call log_warning('The ZBL is true. In the current MiLaDy implementation the 2B kernel should be activated.'//nwl//&
               '... and is not the case.  The master of puppets have decided for you to switch on the  2B kernels by activate_k2b = .true.')
            activate_k2b = .true.

         end if
      end if

      if (zbl_type == zbl_mode_alone) then
         ! zbl_mode_alone: additive mode E_tot = E_ZBL + E_MB
         ! The many-body descriptor uses r_cut_in inner cutoff so E_MB = 0 for r < r_cut_in.
         ! No k2b bridge needed. k2b can still be used independently if the user wants.
         if (r_cut_in > 0.d0) then
            call log_info('ZBL type 2 (additive): E_tot = E_ZBL + E_MB. The many-body part is zero for r < r_cut_in = '//vtoa(r_cut_in))
         else
            call log_info('ZBL type 2 (additive): E_tot = E_ZBL + E_MB. Pair-specific r_cut_in from covalent radii.')
         end if
      end if

   end if  !zbl_potential

   if (activate_k2b) then
      if (r_cut_2b <=0 ) then
         r_cut_2b = r_cut
         r_cut_width_2b = r_cut_width_2b
      end if
   end if


   if (optimize_weights_db) then
      if (lambda_krr > 0.d0) then
         lambda_krr_fake = lambda_krr
      end if

      if (lambda_krr < 0.d0) then
         lambda_krr = (min_lambda_krr + max_lambda_krr)/2.d0
         call log_info('The initial lambda_krr is negative and set-up to the average value in [lambda_krr_min, lambda_krr_max]'//vtoa(lambda_krr))
      end if
      if ((min_lambda_krr < 0.d0) .or. (max_lambda_krr < 0.d0)) then
         call mld_mpi_abort('read_ml: min_lambda_krr or max_lambda_krr < 0 '//vtoa(min_lambda_krr)//vtoa(max_lambda_krr))
      end if
      if (optimize_weights_L1 .and. optimize_weights_L2) then
         call mld_mpi_abort('read_ml: optimize_weigths_L1 and optimize_weigths_L2 cannot be true simulanously '//vtoa(optimize_weights_L1)//vtoa(optimize_weights_L2))
      end if
      if (optimize_weights_L1 .and. optimize_weights_Le) then
         call mld_mpi_abort('read_ml: optimize_weigths_L1 and optimize_weigths_Le cannot be true simulanously '//vtoa(optimize_weights_L1)//vtoa(optimize_weights_Le))
      end if
      if (optimize_weights_L2 .and. optimize_weights_Le) then
         call mld_mpi_abort('read_ml: optimize_weigths_L2 and optimize_weigths_Le cannot be true simulanously '//vtoa(optimize_weights_L2)//vtoa(optimize_weights_Le))
      end if
      if (min_lambda_krr > max_lambda_krr) then
         call mld_mpi_abort('read_ml: lambda_krr_min should be <= lambda_krr_max'//vtoa(min_lambda_krr)//vtoa(max_lambda_krr))
      end if
   end if

!------------------------------------
   if (weighted_3ch) then
      if (.not. (weighted)) then
         call log_warning('read_ml: weighted cannot be false as long as weighted_3ch is true'//nwl// &
            'fixed as set weighted True')
      end if
      weighted = .true.
   end if


   if  ((img_weighted).and.(weighted)) then
      call log_warning('read_ml: weighted cannot be true in the same time as img_weighted '//nwl// &
         'on of them should be fixed to .false. ')
      call mld_mpi_abort('read_ml: milady will stop weighted and img_weighted true in the same time')
      !call mld_mpi_abort('read_ml: milady will stop')
   end if

   if (descriptor_type /= descriptor_milady) then
      if (img_weighted) then
         call log_warning('read_ml: img_weighted is true for other descriptor than 77 '//nwl// &
            'we will put img_weighted = .false. ')
         img_weighted = .false.
      end if
   end if
!--------------------------------------------------------------------------------

   call periodic_table
   call fix_type_of_atoms

   ! Setup pair-specific cutoff arrays (after species are known)
   call setup_pair_cutoffs()

   if ((ml_type == ml_type_krr) .or. (ml_type == ml_type_basis)) then
      if (write_kernel_matrix) then
         call log_warning('read_ml: ml_type is 1 and write_kernel_matrix is True.'//nwl// &
            'kernel cannot be written, fixed as set write_kernel_matrix False')
         write_kernel_matrix = .false.
      end if

      if (ml_type == ml_type_krr) then
         mld_order = mld_kernel
      end if
   end if

   if (ml_type == ml_type_nlinear) then
      mld_order = mld_linear_extended
   end if



   select case (mld_order)

    case (mld_linear)
      call log_info('Linear ML')
    case (mld_linear_extended)
      call log_info('n-Linear ML')
    case (mld_kernel)
      call log_info('Kernel + linear precondition ML')
    case (mld_quadratic)
      call log_info('Quadratic ML')
      select case (mld_type_quadratic)
       case (mld_type_quadratic_qnml)
         call log_info('QML - LINEAR PRECONDITION')
       case (mld_type_quadratic_qml)
         call log_info('QML - FULL')
       case (mld_type_quadratic_bilinear)
         call log_info('QML - BI-LINEAR')
       case (mld_type_quadratic_zaxa)
          call log_info('QML - zaxa')
       case (mld_type_quadratic_ZX)
          call log_info('QML - ZX')   
       case default
         call mld_mpi_abort('mld_type_quadratic only 1, 2, 3, 5 and 6. Read the manual.')
      end select

    case (mld_polyc)
      call log_info('Poly Chaos ML')
      if (polyc_n_hermite < 1) then
         call mld_mpi_abort('read_ml: Poly_Chaos ML requires Hermite polynomials order >= 2, polyc_n_hermite: '//vtoa(polyc_n_hermite)//nwl// &
            'fix as change polyc_n_hermite in ml file')
      end if
      if (polyc_n_hermite > 4) then
         call mld_mpi_abort('read_ml: Poly_Chaos ML requires Hermite polynomials order <= 4, polyc_n_hermite: '//vtoa(polyc_n_hermite)//nwl// &
            'fix as change polyc_n_hermite in ml file')
      end if
      if (polyc_n_poly > 3) then
         call mld_mpi_abort('read_ml: Poly_Chaos ML requires the polynomial order <= 3, polyc_n_poly: '//vtoa(polyc_n_poly)//nwl// &
            'fix as change polyc_n_poly in ml file')
      end if
   end select

   if (.not. ((mld_order == mld_linear) .or. (mld_order == mld_quadratic) &
      .or. (mld_order == mld_polyc) &
      .or. (mld_order == mld_kernel) &
      .or. (mld_order == mld_linear_extended))) then
      call mld_mpi_abort('read_ml: snap potential out of range, mld_order: '//vtoa(mld_order))
   end if

   !TODO: to_remove  for old world .... 
   nd_data=1


   if (.not. ((mld_fit_type == 0) .or. &
      (mld_fit_type == 1) .or. &
      (mld_fit_type == 2) .or. &
      (mld_fit_type == 10) .or. &
      (mld_fit_type == 4) .or. &
      (mld_fit_type == 5) .or. &
      (mld_fit_type == 6) .or. &
      (mld_fit_type == 3))) then

      call mld_mpi_abort('read_ml: type of fitting out of range, mld_fit_type: '//vtoa(mld_fit_type)//nwl// &
         'should be 0, 1, 2, 3, 4, 5, 6 or 10')
   end if

   if (mld_fit_type == 2 ) then 
     call log_critical("fit with constraints is not available in this MLD version")
     call mld_mpi_abort('read_ml: abort for this mld_fit_type  2. Change it ! ')
   end if 


   if (kcross .and. marginal_likelihood) call mld_mpi_abort('read_ml: both kcross or marginal_likelihood cannot be true')

   if (ml_type /= ml_type_descriptors) then
      if (.not. (desc_forces)) call log_warning('read_ml: forces are desactivated. if you want forces set desc_forces True')
   end if


   if (ml_type == ml_type_descriptors) then
      train_only = .true.
      if (.not.(write_desc)) then
         write_desc=.true.
         call log_warning('read_ml: write_desc =.false.. In this ML mode writing descriptors are activated even if they were originally .false.')
      end if
   end if

   if (ml_type /= ml_type_descriptors) then
      if (lmask) then
         call log_warning('read_ml: lmask =.true. In this ML mode lmask can be only .false. So will be set to false')
         lmask = .false.
      end if

      if (mask_file) then
         call log_warning('read_ml: mask_file =.true. In this ML mode mask_file can be only .false. So will be set to false')
         mask_file = .false.

      end if
   end if

   if (read_desc_dump .or. write_desc_dump) then
      if (desc_file_format /= 1 ) then
         call log_warning('read_ml: desc_file_format /=1 in those read_desc_dump and write_desc_dump mode')
         call log_warning('read_ml: can be only 1. it wil be switched to 1.')
         desc_file_format = 1
      end if
   end if


   if ((descriptor_type == descriptor_behler) .or. &
      (descriptor_type == descriptor_g2) .or. &
      (descriptor_type == descriptor_g3) .or. &
      (descriptor_type == descriptor_g2_pow_so4) .or. &
      (descriptor_type == descriptor_g2_afs) .or. &
      (descriptor_type == descriptor_g2_bispectrum_so4)) then
      if (strict_behler) then
         if (descriptor_type /= descriptor_behler) then
            call log_warning('read_ml: if strict_behler is True, descriptor should be Behler as descriptor_type 3')
         end if
      end if
   end if

   if ((descriptor_type == descriptor_g2) .and. ((n_g2_eta <= 0) .or. (n_g2_rs <= 0))) then
      call log_warning('n_g2_eta and n_g2_rs should be > 0')
   elseif ((descriptor_type == descriptor_g3) .and. ((n_g3_eta <= 0) .or. (n_g3_zeta <= 0) .or. (n_g3_lambda <= 0))) then
      call log_warning('n_g3_eta,n_g3_zeta,n_g3_lambda should be > 0')
   elseif ((descriptor_type == descriptor_behler) .and. ((n_g2_eta <= 0) .or. (n_g2_rs <= 0) .or. (n_g3_zeta <= 0) .or. (n_g3_lambda <= 0))) then
      call log_warning('n_g2 g3_eta,n_g2_rs,n_g3_zeta,n_g3_lambda should be > 0')
   elseif ((descriptor_type == descriptor_afs) .and. ((n_rbf <= 0) .or. (n_cheb < 0))) then
      call log_warning('n_rbf should be > 0 and n_cheb should be >= 0')
   elseif ((descriptor_type == 6) .and. ((n_rbf <= 0) .or. (l_max < 0))) then
      call log_warning('n_rbf should be > 0 and l_max should be >= 0')
   elseif ((descriptor_type == 7) .and. ((n_rbf <= 0) .or. (l_max < 0))) then
      call log_warning('n_rbf should be > 0 and l_max should be >= 0')
   end if

   
   if (descriptor_type == descriptor_ace) then
       if (.not.((ace_chem == 0) .or. (ace_chem == 1) .or. (ace_chem == 2))) then
          call log_critical('read_ml: Consistency check failed')
          call log_critical('ace_chem should be 0 1 or 2 , but is '//vtoa(ace_chem))
          call mld_mpi_abort('read_ml: ace_chem should be 0, 1 or 2')
       end if

       if (ace_chem == 2 ) then
           ace_radial_chem = 2 
       end if 
       ntest = nitems2(ace_nmax_list)
       if (ntest .lt. ace_numax) then
          call log_critical('read_ml: Consistency check failed')
          call log_critical('#members of ace_nmax_list /= nu_ace_max '//vtoa(ntest)//vtoa(ace_numax))
          call log_critical("ace_nmax_list should have the same number of terms as the value of ace_numax")
          call log_critical("ace_nmax_list = "//" "//ace_nmax_list)
          call log_critical("nu_ace_max = "//vtoa(ace_numax))
          call mld_mpi_abort('read_ml: #members of ace_nmax_list  /= ace_numax')
       end if
       if (allocated(ace_nmax)) deallocate(ace_nmax) ; allocate(ace_nmax(ace_numax))
       read(ace_nmax_list,*) ace_nmax

       ntest = nitems2(ace_kmax_list)
       if (ntest .lt. ace_numax) then
          call log_critical('read_ml: Consistency check failed')
          call log_critical('#members of ace_kmax_list /= nu_ace_max '//vtoa(ntest)//vtoa(ace_numax))
          call log_critical("ace_kmax_list should have the same number of terms as the value of ace_numax")
          call log_critical("ace_kmax_list = "//" "//ace_kmax_list)
          call log_critical("nu_ace_max = "//vtoa(ace_numax))
          call mld_mpi_abort('read_ml: #members of ace_kmax_list  /= ace_numax')
       end if
       if (allocated(ace_kmax)) deallocate(ace_kmax) ; allocate(ace_kmax(ace_numax))
       read(ace_kmax_list,*) ace_kmax

       ntest = nitems2(ace_lmax_list)
       if (ntest .lt.  ace_numax) then
          call log_critical('read_ml: Consistency check failed')
          call log_critical('#members of ace_lmax_list /= nu_ace_max '//vtoa(ntest)//vtoa(ace_numax))
          call log_critical("ace_lmax_list should have the same number of terms as the value of ace_numax")
          call log_critical("ace_lmax_list = "//" "//ace_lmax_list)
          call log_critical("nu_ace_max = "//vtoa(ace_numax))
          call mld_mpi_abort('read_ml: #members of ace_lmax_list  /= ace_numax')
       end if
       if (allocated(ace_lmax)) deallocate(ace_lmax) ; allocate(ace_lmax(ace_numax))
       read(ace_lmax_list,*) ace_lmax

       ntest = nitems2(ace_lambda_list)
       if (ntest .lt. ace_numax) then
          call log_critical('read_ml: Consistency check failed')
          call log_critical('#members of ace_lambda_list /= nu_ace_max '//vtoa(ntest)//vtoa(ace_numax))
          call log_critical("ace_lambda_list should have the same number of terms as the value of ace_numax")
          call log_critical("ace_lambda_list = "//" "//ace_lambda_list)
          call log_critical("nu_ace_max = "//vtoa(ace_numax))
          call mld_mpi_abort('read_ml: #members of ace_lambda_list  /= ace_numax')
       end if
       if (allocated(ace_lambda)) deallocate(ace_lambda) ; allocate(ace_lambda(ace_numax))
       read(ace_lambda_list,*) ace_lambda

       if (zetaace_order .gt. 1) then 
         if (zetaace_order .gt. MAX_ZETAACE_ORDER) then 
            call log_critical('read_ml: Consistency check failed')
            call log_critical('zetaace_order > MAX_ZETAACE_ORDER '//vtoa(zetaace_order)//vtoa(MAX_ZETAACE_ORDER))
            call log_critical('zetaace_order should be <= MAX_ZETAACE_ORDER = '//vtoa(MAX_ZETAACE_ORDER))
            call mld_mpi_abort('read_ml: zetaace_order > MAX_ZETAACE_ORDER')
         end if
       end if 

   end if

   if (descriptor_type == descriptor_zetabody) then 
      if (l_body_order(1)) then 
         call log_warning('l_body_order(1) is True. For ZETA_BODY should be false. We will turn to False')
         l_body_order(1) = .false. 
      end if 

      if (l_body_order(2)) then 
        call log_warning('l_body_order(2) is true. For ZETA_BODY, it is recommended to set this to False.')
        call log_warning('It is possible to use 2-body interaction, but it is recommended to use active_k2b =.true. instead.') 
        !$! l_body_order(2) = .false. 
      end if 

      if (l_body_order(4)) then 
        call log_warning('l_body_order(4) is true. For ZETA_BODY this should be set to False. We will do it.')
        l_body_order(4) = .false. 
      end if 

      if (l_body_order(5)) then 
        call log_warning('l_body_order(5) is true. For ZETA_BODY this should be set to False. We will do it.')
        l_body_order(5) = .false. 
      end if 

      if (l_body_order(6)) then 
        call log_warning('l_body_order(6) is true. For ZETA_BODY this should be set to False. We will do it.')
        l_body_order(6) = .false. 
      end if 

      if (r_cut_z2b .lt. 0) r_cut_z2b = r_cut 
      if (r_cut_z3b .lt. 0) r_cut_z3b = r_cut

      if (r_cut_width_z2b .lt. 0) r_cut_width_z2b = r_cut_width  
      if (r_cut_width_z3b .lt. 0) r_cut_width_z3b = r_cut_width  

      if (zetabody_order .ge. MAX_ZETABODY_ORDER) then 
         call log_critical('read_ml: Consistency check failed')
         call log_critical('zetabody_order >= MAX_ZETABODY_ORDER '//vtoa(zetabody_order)//vtoa(MAX_ZETABODY_ORDER))
         call log_critical('zetabody_order should be < MAX_ZETABODY_ORDER = '//vtoa(MAX_ZETABODY_ORDER))
         call mld_mpi_abort('read_ml: zetabody_order >= MAX_ZETABODY_ORDER')
      end if

   end if 

   if (descriptor_type == descriptor_ftnbody) then
      !TODO check_dimensions ....
      read (dim_fourier_nbody, *) dim_rff(:)
      read (length_fourier_nbody, *) length_rff(:)
      read (delta_fourier_nbody, *) delta_rff(:)
      ! multispecies chemistry-mode parameters (docs/perspective_ftnbody.md)
      read (ftnbody_chem_rank, *) ftnbody_chem_rank_n(:)
      read (ftnbody_hash_channels, *) ftnbody_hash_channels_n(:)
      ! sPIP maximum polynomial degree per body order (docs/MLT5 sec. 2.2.6)
      read (spip_degree_nbody, *) spip_degree_n(:)
      ! geometry model (docs/MLT5 sec. 2.2.4 and 2.2.6):
      !   'poly'/'pip' compact invariant q (default), 'gramm' ordered Gram
      !   coordinates + orbit averaging, 'cpip' compact mixed invariant set,
      !   'spip' systematic orbit-polynomial set
      select case (trim(adjustl(ftnbody_model)))
      case ('poly', 'POLY', 'Poly', 'pip', 'PIP', 'Pip')
         ftnbody_model_id = FTNBODY_MODEL_POLY
      case ('gramm', 'GRAMM', 'Gramm', 'gram', 'GRAM', 'Gram')
         ftnbody_model_id = FTNBODY_MODEL_GRAMM
      case ('cpip', 'CPIP', 'Cpip', 'cPIP')
         ftnbody_model_id = FTNBODY_MODEL_CPIP
      case ('spip', 'SPIP', 'Spip', 'sPIP')
         ftnbody_model_id = FTNBODY_MODEL_SPIP
      case default
         call mld_critical_abort('read_ml: unknown ftnbody_model "'//trim(adjustl(ftnbody_model))// &
                                 '" (allowed: poly/pip, gramm, cpip, spip)')
      end select
      ! cPIP is deliberately restricted to the chemistry-blind and low-rank modes
      if (ftnbody_model_id == FTNBODY_MODEL_CPIP .and. &
          ftnbody_chem_mode /= 0 .and. ftnbody_chem_mode /= 2) then
         call mld_critical_abort('read_ml: ftnbody_model "cpip" supports only ftnbody_chem_mode 0 or 2')
      end if
      call log_info('ML: FT-nBody geometry model: '//trim(adjustl(ftnbody_model)))
      if ((l_body_order(1)) .and. (dim_rff(1) == 0)) then
         call log_warning('l_body_oder(1) is True, yet dim_rff(1) is 0, we set default value to np_kernel_full :'//vtoa(np_kernel_full))
         if (np_kernel_full == 0) then
            call mld_critical_abort('read_ml: Consistency check failed' //nwl// &
               'l_body_oder(1) is True, yet dim_rff(1) is 0 and np_kernel_full is 0')
         else
            dim_rff(1) = np_kernel_full
         end if
      end if
      if ((l_body_order(2)) .and. (dim_rff(2) == 0)) then
         call log_warning('l_body_oder(2) is True, yet dim_rff(2) is 0, we set default value to np_kernel_full :'//vtoa(np_kernel_full))
         if (np_kernel_full == 0) then
            call mld_critical_abort('read_ml: Consistency check failed' //nwl// &
               'l_body_oder(2) is True, yet dim_rff(2) is 0 and np_kernel_full is 0')
         else
            dim_rff(2) = np_kernel_full
         end if
      end if
      if ((l_body_order(3)) .and. (dim_rff(3) == 0)) then
         call log_warning('l_body_oder(3) is True, yet dim_rff(3) is 0, we set default value to np_kernel_full :'//vtoa(np_kernel_full))
         if (np_kernel_full == 0) then
            call mld_critical_abort('read_ml: Consistency check failed' //nwl// &
               'l_body_oder(3) is True, yet dim_rff(3) is 0 and np_kernel_full is 0')
         else
            dim_rff(3) = np_kernel_full
         end if
      end if
      if ((l_body_order(4)) .and. (dim_rff(4) == 0)) then
         call log_warning('l_body_oder(4) is True, yet dim_rff(4) is 0, we set default value to np_kernel_full :'//vtoa(np_kernel_full))
         if (np_kernel_full == 0) then
            call mld_critical_abort('read_ml: Consistency check failed' //nwl// &
               'l_body_oder(4) is True, yet dim_rff(4) is 0 and np_kernel_full is 0')
         else
            dim_rff(4) = np_kernel_full
         end if
      end if
      if ((l_body_order(5)) .and. (dim_rff(5) == 0)) then
         call log_warning('l_body_oder(5) is True, yet dim_rff(5) is 0, we set default value to np_kernel_full :'//vtoa(np_kernel_full))
         if (np_kernel_full == 0) then
            call mld_critical_abort('read_ml: Consistency check failed' //nwl// &
               'l_body_oder(5) is True, yet dim_rff(5) is 0 and np_kernel_full is 0')
         else
            dim_rff(5) = np_kernel_full
         end if
      end if
   end if


   if (selection_type == selection_type_first ) then
      call log_info('ML: selection of the first no selected configs of the database class / db line ')
   elseif (selection_type == selection_type_last) then
      call log_info('ML: selection of the last no selected configs of the database class / db line')
   elseif (selection_type == selection_type_random) then
      call log_info('ML: selection of random no selected configs of the database class / db line ')
   elseif (selection_type == selection_type_first_start) then
      call log_info('ML: selection starts not from 1 but from a given  start in db line ')   
   else
      call mld_mpi_abort('read_ml: selection_type should be 1,2 or 3: '//vtoa(selection_type))
   end if

   if (.not. (desc_forces)) then
      call log_warning('descriptors for forces are not computed as desc_forces False.'//nwl// &
         'you do WHATEVER YOU WANT ... even without THE FORCE')
   end if

   if (train_only) then
      call log_warning('only TRAINING configuration. If you want to test also set train_only False')
   end if

   if (kcross) then
      if (n_kcross <= 0) then
         call mld_mpi_abort('read_ml: n_kcross should > 0. Good choice is in [5, 10]'//nwl// &
            'but it can be any value. You figure out that the training time will be x n_kcross')
      end if
      if (nd_data <= 0) then
         call mld_mpi_abort('read_ml: n_data (size if the database) should be > 0')
      end if
   end if

   select case (ml_type)
    case (0)
      call log_info('Machine Learning type with linear basis, LML')
    case (1)
      call log_info('Machine learning type with Kernel Ridge Regression')
    case (2)
      call mld_mpi_abort('read_ml: Machine learning type with Gaussian Precesses (not implemented, yet)')
   end select

   if (write_design_matrix) then
      if (ml_type < 0) then
         call mld_mpi_abort('read_ml: for dump if write_design_matrix True, ml_type should be >= 0')
      end if
   end if

   if (write_test_design_matrix) then
      if (ml_type < 0) then
         call mld_mpi_abort('read_ml: for dump if write_test_design_matrix True, ml_type should be >= 0')
      end if
   end if


   if (write_kernel_matrix) then
      if (ml_type /= ml_type_analysis) then
         call log_warning('read_ml: if write_kernel_matrix=.true. and dump of the kernel can be done only in the model ml_type = -2')
         call log_warning('read_ml: was write_kernel_matrix=.true. and now it bacomes write_kernel_matrix=.false. ')
         write_kernel_matrix = .false.
      end if
   end if


   if ((ml_type .eq. 1) .or. (ml_type .eq. 2)) then
      select case (kernel_type)
       case (kernel_se)
         call log_info('kernel SQUARED EXPONENTIAL sigma_kernel '//vtoa(sigma_kernel)//'length_kernel '//vtoa(length_kernel))
       case (kernel_ou)
         call mld_mpi_abort('kernel ORSTEIN UHLENBECK not yet implemented')
       case (kernel_mc)
         call mld_mpi_abort('kernel MATERN CLASS not yet implemented')
       case (kernel_po)
         call log_info('kernel POLYNOMIAL sigma_kernel '//vtoa(sigma_kernel)//'length_kernel '//vtoa(length_kernel)//'Poly order kernel '//vtoa(kernel_power))
       case (kernel_po_scaled)
         call log_info('kernel POLYNOMIAL SCALED sigma_kernel '//vtoa(sigma_kernel)//'Poly order kernel '//vtoa(kernel_power))
       case (kernel_maha)
         call log_info('kernel MAHA MILADY')
       case (kernel_random_maha)
         call log_info('kernel RANDOM MAHA MILADY')
       case (kernel_random)
         call log_info('kernel RANDOM Gaussian MILADY')
       case (kernel_random_po)
         call log_info('kernel RANDOM POLY MILADY')
       case default
         call mld_mpi_abort('No implementation for kernel_type '//vtoa(kernel_type))
      end select
   end if

   if (toy_model) then
      call log_info('TOY model - ML_perfect descriptors')
   else

      append_infos = ""
      if (activate_k2b) then
         append_infos = trim(append_infos)//" kernel 2-body +"
      end if

      select case (descriptor_type)
       case (descriptor_g2)
         call log_info('NDM + '//trim(append_infos)//' ML_G2')
       case (descriptor_g3)
         call log_info('NDM + '//trim(append_infos)//' ML_G3')
       case (descriptor_behler)
         call log_info('NDM + '//trim(append_infos)//' ML_BEHLER')
       case (descriptor_afs)
         call log_info('NDM + '//trim(append_infos)//' ML_AFS')
       case (descriptor_g2_afs)
         call log_info('NDM + '//trim(append_infos)//' ML_G2_AFS')
       case (descriptor_pow_so3)
         call log_info('NDM + '//trim(append_infos)//' ML_POW_SO3')
       case (descriptor_pow_so3_3body)
         call log_info('NDM + '//trim(append_infos)//' ML_POW_SO3 3BODY')
       case (descriptor_pow_so4)
         call log_info('NDM + '//trim(append_infos)//' ML_POW_SO4')
       case (descriptor_g2_pow_so4)
         call log_info('NDM + '//trim(append_infos)//' ML_G2_POW_SO4')
       case (descriptor_bispectrum_so4)
         call log_info('NDM + '//trim(append_infos)//' ML_BISPECTRUM_SO4')
       case (descriptor_mtp)
         call log_info('NDM + '//trim(append_infos)//' MTP')
       case (descriptor_g2_bispectrum_so4)
         call log_info('NDM + '//trim(append_infos)//' ML_G2_BISPECTRUM_SO4')
       case (descriptor_milady)
         call log_info('NDM + '//trim(append_infos)//' ML_MILADY')
       case (descriptor_body)
         call log_info('NDM + '//trim(append_infos)//' ML_BODY single element')
       case (descriptor_zetabody)
         call log_info('NDM + '//trim(append_infos)//' ML_ZETA_BODY direct space')
       case (descriptor_ace)
         call log_info('NDM + '//trim(append_infos)//' standard ACE')
       case (descriptor_ftnbody)
         call log_info('NDM + '//trim(append_infos)//' ML_FOURIER_NBODY N BODY FOURIER SPACE')
       case (descriptor_tbind)
         call log_info('NDM + '//trim(append_infos)//' ML_TBIND_NBODY N BODY TIGHT BINDING')
       case default
         call mld_mpi_abort('message from read_ml_file: No implementation for descriptor_type '//vtoa(descriptor_type))
      end select
   end if

   if (optimize_weights_db) then
      none_class = nitems(class_no_optimize_weights)
      call log_info('Weights optimization will be performed')
      if (allocated(no_class_weights)) deallocate (no_class_weights); allocate (no_class_weights(none_class))
      read (class_no_optimize_weights, *) no_class_weights
      tmp = 'no weights optimization for classes:'//nwl
      do i = 1, none_class
         tmp = tmp//vtoa(i)//no_class_weights(i)//nwl
      end do
      call log_info(tmp)
   end if

   if (optimize_weights_db .eqv. optimize_weights_chem) then
      if (optimize_weights_chem .eqv. .true.) then
         call mld_mpi_abort('FATAL: optimize_weights_db and optimize_weights_chem cannot be true in the same time. At least one should be .false.')
      end if
   end if

   if (ml_type == ml_type_analysis) then

      call log_info('-- MCD / Kernel parameters --')
      int_local = nitems2(classes_for_mcd)
      if (int_local == 0) then
         call mld_mpi_abort('the list of classes_for_mcd is empty. Please put at least one class'//nwl// &
            'an example classes_for_mcd = "01 05 O6"')
      end if
      call log_info('number of classes for MCD '//vtoa(int_local)//nwl// &
         'classes for MCD '//classes_for_mcd)
      if (allocated(classes_train_for_sigma)) deallocate (classes_train_for_sigma)
      allocate (classes_train_for_sigma(int_local))
      read (classes_for_mcd, *) classes_train_for_sigma(:)



   end if



   if ((write_desc .or. write_desc_dump) .and. read_desc_dump) then
      if (rangml == 0) then
         tmp = 'Confusing choice in read/write parameters.'//nwl// &
            'You should decide if this run is for reading or writing'//nwl// &
            'Currently you have choose:'//nwl

         tmp = tmp//toline('write_desc', vtoa(write_desc))//nwl
         tmp = tmp//toline('write_desc_dump', vtoa(write_desc_dump))//nwl
         tmp = tmp//toline('read_desc_dump', vtoa(read_desc_dump))//nwl

         call log_info(tmp)
      end if
      text = "Confuse read and write option in read_ml"

      !call end_ml(text, scalapack_driver, context, rangml)
   end if

   call write_all_input

end subroutine read_ml_file


subroutine set_grid_lambda_krr()

   use ml_in_ndm_module, only: min_lambda_krr, max_lambda_krr, vector_lambda_krr, n_values_lambda_krr, lambda_krr
   use mesh_grid, only: log_grid
   use mld_string
   use mld_logger
   use mld_mpi

   implicit none

   if (n_values_lambda_krr .lt. 0) then
      call mld_mpi_abort('n_values_lambda_krr cannot be negative: '//vtoa(n_values_lambda_krr))
   end if

   if (n_values_lambda_krr == 1) then
      if ((min_lambda_krr /= max_lambda_krr)) then
         call mld_mpi_abort('if n_values_lambda_krr = 1, requires min_lambda_krr = max_lambda_krr: '// &
            vtoa(min_lambda_krr)//vtoa(max_lambda_krr))
      end if
   end if



   if ((min_lambda_krr .lt. 0) .or. (max_lambda_krr .lt. 0)) then
      n_values_lambda_krr = 21
      if (allocated(vector_lambda_krr)) deallocate (vector_lambda_krr)
      allocate (vector_lambda_krr(n_values_lambda_krr))
      vector_lambda_krr = (/1.d-10, 1.d-09, 1.d-08, 1.d-07, &
         1.d-06, 1.d-05, 1.d-04, 1.d-03, &
         1.d-02, 1.d-01, 1.d+00, 1.d+01, &
         1.d+02, 1.d+03, 1.d+04, 1.d+05, &
         1.d+06, 1.d+07, 1.d+08, 1.d+09, &
         1.d+10/)
   else if ((min_lambda_krr == max_lambda_krr)) then
      n_values_lambda_krr = 1
      if (allocated(vector_lambda_krr)) deallocate (vector_lambda_krr)
      allocate (vector_lambda_krr(n_values_lambda_krr))
      vector_lambda_krr = (/min_lambda_krr/)
      lambda_krr = min_lambda_krr
   else
      !if (allocated(vector_lambda_krr)) deallocate (vector_lambda_krr)
      !allocate (vector_lambda_krr(n_values_lambda_krr))
      call log_grid(n_values_lambda_krr, min_lambda_krr, max_lambda_krr, vector_lambda_krr)
      !delta_lambda = (log10(max_lambda_krr) - log10(min_lambda_krr))/dble(n_values_lambda_krr - 1)
      !do i = 1, n_values_lambda_krr
      !  vector_lambda_krr(i) = 10.d0**(log10(min_lambda_krr) + dble(i - 1)*delta_lambda)
      !end do
   end if

end subroutine set_grid_lambda_krr


subroutine fix_type_of_atoms

   use ml_in_ndm_module, only: fix_weighted_for_element, fix_weighted_for_element_3ch, &
      fix_weighted_for_element_ini, fix_weighted_for_element_3ch_ini, &
      weighted_auto, &
      fix_no_of_elements_invisible, fix_ch_elements_invisible, fix_mass_elements_invisible, &
      fix_type_to_periodic_invisible, fix_Z_elements_invisible, &
      weighted, weighted_3ch, chemical_elements_invisible, linvisible

   use module_chemical_species, only: periodic_table_element, size_periodic_table, &
      map_species_full, map_species_half, &
      size_species_full, size_species_half, &
      fix_type_to_periodic, fix_covalent_radius_elements, fix_mass_elements, &
      fix_Z_elements, fix_ch_elements, fix_no_of_elements, &
      chemical_elements,  weight_per_element, weight_per_element_3ch, &
      img_num_ch, img_weighted,  ws1, ws2, ws3, ws4, fix_wspecies

   use math
   use module_db_poscar, only: fix_ref_energy_per_element, ref_energy_per_element
   use derived_types, only: typ_species_full, typ_species_half
   use mld_mpi
   use mld_string
   use mld_logger

   implicit none

   character(len=:), allocatable    :: tmp
   integer  :: int_local, nitems2, i_p, icount
   character(len=1)     :: quote, dquote
   integer :: ii, jj

   quote = char(39)
   dquote = char(34)

   int_local = nitems2(chemical_elements)
   if (int_local == 0) then
      call mld_mpi_abort("fix_type_of_atoms: read_ml_file with chemical_elements"//nwl// &
         "the list of chemical elements is empty"//nwl// &
         "please put the appropiate chemical_elements list"//nwl// &
         'example: chemical_elements="Fe W Os"')
   end if


   if (int_local /= fix_no_of_elements) then
      call mld_mpi_abort('the list of chemical elements and the list of weigths for each elements has not the same lengths'//nwl// &
         'case A expected lenght: '//vtoa(fix_no_of_elements))
   end if


   fix_no_of_elements = int_local
   if (allocated(fix_ch_elements)) deallocate (fix_ch_elements)
   allocate (fix_ch_elements(int_local))
   if (allocated(fix_Z_elements)) deallocate (fix_Z_elements)
   allocate (fix_Z_elements(int_local))
   if (allocated(fix_mass_elements)) deallocate (fix_mass_elements)
   allocate (fix_mass_elements(int_local))
   if (allocated(fix_covalent_radius_elements)) deallocate (fix_covalent_radius_elements)
   allocate (fix_covalent_radius_elements(int_local))
   if (allocated(fix_type_to_periodic)) deallocate (fix_type_to_periodic)
   allocate (fix_type_to_periodic(int_local))
   read (chemical_elements, *) fix_ch_elements(:)

   icount = 0
   tmp = "-- DATABASE ATOMS --"//nwl// &
      " id element   Z     mass"//nwl
   do int_local = 1, size(fix_ch_elements, dim=1)
      do i_p = 1, size_periodic_table
         if (periodic_table_element(i_p)%symbol == fix_ch_elements(int_local)) then
            icount = icount + 1
            fix_Z_elements(int_local) = periodic_table_element(i_p)%Z
            fix_mass_elements(int_local) = periodic_table_element(i_p)%mass
            fix_type_to_periodic(int_local) = i_p
            fix_covalent_radius_elements(int_local) = periodic_table_element(i_p)%covalent_radius
         end if
      end do
      tmp = tmp//vtoa(int_local, 3)//vtoa(fix_type_to_periodic(int_local), 4)//fix_ch_elements(int_local)// &
         vtoa(fix_Z_elements(int_local))//vtoa(fix_mass_elements(int_local))//nwl
   end do
   call log_info(tmp)

   ! if (ml_type >= 0 ) then
   int_local = nitems2(ref_energy_per_element)
   if (int_local /= fix_no_of_elements) then
      call log_warning('the list of reference energy and the list of no of elements has not the same lengths'//nwl// &
         'current and expected are '//vtoa([int_local, fix_no_of_elements])//nwl// &
         'fixed for you as ALL ref_energy will be zero')


      if (allocated(fix_ref_energy_per_element)) deallocate (fix_ref_energy_per_element)
      allocate (fix_ref_energy_per_element(fix_no_of_elements))
      fix_ref_energy_per_element(:) = 0.d0
      ! call mld_mpi_abort('problem length ref_energy_per_element"
   else
      if (allocated(fix_ref_energy_per_element)) deallocate (fix_ref_energy_per_element)
      allocate (fix_ref_energy_per_element(fix_no_of_elements))
      read (ref_energy_per_element, *) fix_ref_energy_per_element(:)
   end if
   ! end if


   if (img_weighted) then
      if (allocated(fix_wspecies)) deallocate(fix_wspecies) ; allocate(fix_wspecies(fix_no_of_elements,img_num_ch))
      select case (img_num_ch)

       case(1)
         int_local = nitems2(ws1)
         if (int_local /= fix_no_of_elements) then
            call mld_mpi_abort('list of chemical elements and the list of weigths ws1 for each elements has not the same lengths'//nwl// &
               'case B expected lenght: '//vtoa(fix_no_of_elements))
         end if
         read (ws1, *) fix_wspecies(:,1)

       case (2)
         int_local = nitems2(ws1)
         if (int_local /= fix_no_of_elements) then
            call mld_mpi_abort('list of chemical elements and the list of weigths ws1 for each elements has not the same lengths'//nwl// &
               'case B expected lenght: '//vtoa(fix_no_of_elements))
         end if
         int_local = nitems2(ws2)
         if (int_local /= fix_no_of_elements) then
            call mld_mpi_abort('list of chemical elements and the list of weigths ws2 for each elements has not the same lengths'//nwl// &
               'case B expected lenght: '//vtoa(fix_no_of_elements))
         end if
         read (ws1, *) fix_wspecies(:,1)
         read (ws2, *) fix_wspecies(:,2)


       case(3)
         int_local = nitems2(ws1)
         if (int_local /= fix_no_of_elements) then
            call mld_mpi_abort('list of chemical elements and the list of weigths ws1 for each elements has not the same lengths'//nwl// &
               'case B expected lenght: '//vtoa(fix_no_of_elements))
         end if
         int_local = nitems2(ws2)
         if (int_local /= fix_no_of_elements) then
            call mld_mpi_abort('list of chemical elements and the list of weigths ws2 for each elements has not the same lengths'//nwl// &
               'case B expected lenght: '//vtoa(fix_no_of_elements))
         end if
         int_local = nitems2(ws3)
         if (int_local /= fix_no_of_elements) then
            call mld_mpi_abort('list of chemical elements and the list of weigths ws3 for each elements has not the same lengths'//nwl// &
               'case B expected lenght: '//vtoa(fix_no_of_elements))
         end if
         read (ws1, *) fix_wspecies(:,1)
         read (ws2, *) fix_wspecies(:,2)
         read (ws3, *) fix_wspecies(:,3)


       case(4)
         int_local = nitems2(ws1)
         if (int_local /= fix_no_of_elements) then
            call mld_mpi_abort('list of chemical elements and the list of weigths ws1 for each elements has not the same lengths'//nwl// &
               'case B expected lenght: '//vtoa(fix_no_of_elements))
         end if
         int_local = nitems2(ws2)
         if (int_local /= fix_no_of_elements) then
            call mld_mpi_abort('list of chemical elements and the list of weigths ws2 for each elements has not the same lengths'//nwl// &
               'case B expected lenght: '//vtoa(fix_no_of_elements))
         end if
         int_local = nitems2(ws3)
         if (int_local /= fix_no_of_elements) then
            call mld_mpi_abort('list of chemical elements and the list of weigths ws3 for each elements has not the same lengths'//nwl// &
               'case B expected lenght: '//vtoa(fix_no_of_elements))
         end if
         int_local = nitems2(ws4)
         if (int_local /= fix_no_of_elements) then
            call mld_mpi_abort('list of chemical elements and the list of weigths ws4 for each elements has not the same lengths'//nwl// &
               'case B expected lenght: '//vtoa(fix_no_of_elements))
         end if
         read (ws1, *) fix_wspecies(:,1)
         read (ws2, *) fix_wspecies(:,2)
         read (ws3, *) fix_wspecies(:,3)
         read (ws4, *) fix_wspecies(:,4)

       case default
         call mld_mpi_abort('no implementation for this number of channels (max 4 and strict positive): '//vtoa(img_num_ch))

      end select
      weighted = .false.
      weighted_3ch = .false.
      weighted_auto = .false.
   endif

   if (weighted) then
      !
      if (weighted_auto) then
         int_local = fix_no_of_elements
      else
         int_local = nitems2(weight_per_element)
      end if

      if (.not. (weighted_auto)) then
         if (int_local == 0) then
            call mld_mpi_abort('list of weigth_per_element is empty. Please put the appropiate weights')
         else
            if (int_local /= fix_no_of_elements) then
               call mld_mpi_abort('list of chemical elements and the list of weigths for each elements has not the same lengths'//nwl// &
                  'case B expected lenght: '//vtoa(fix_no_of_elements))
            end if
         end if
      end if

      if (allocated(fix_weighted_for_element)) deallocate (fix_weighted_for_element)
      allocate (fix_weighted_for_element(fix_no_of_elements))
      
      if (weighted_auto) then 
       do i_p = 1, fix_no_of_elements
          fix_weighted_for_element(i_p) = fix_Z_elements(i_p)
       end do 
      else 
         read (weight_per_element, *) fix_weighted_for_element(:)
      end if 
      
   end if

   if ((weighted) .and. (weighted_3ch)) then
      !
      if (weighted_auto) then
         int_local = fix_no_of_elements
      else
         int_local = nitems2(weight_per_element_3ch)
      end if


      if (.not. (weighted_auto)) then
         if (int_local == 0) then
            call mld_mpi_abort('list of weigth_per_element_3ch is empty. Please put the appropiate weights')
         else
            if (int_local /= fix_no_of_elements) then
               call mld_mpi_abort('list of chemical elements and the list of weigths 3ch for each elements has not the same lengths'//nwl// &
                  'case D expected lenght: '//vtoa(fix_no_of_elements))
            end if
         end if
      end if

      if (allocated(fix_weighted_for_element_3ch)) deallocate (fix_weighted_for_element_3ch)
      allocate (fix_weighted_for_element_3ch(fix_no_of_elements))
      if (weighted_auto) then 
        do i_p = 1, fix_no_of_elements
          fix_weighted_for_element_3ch(i_p) =  dsqrt(fix_mass_elements(i_p))
          fix_weighted_for_element_3ch(i_p) =  fix_mass_elements(i_p)**2 + fix_covalent_radius_elements(i_p)**2
        end do  
      else 
        read (weight_per_element_3ch, *) fix_weighted_for_element_3ch(:)  
      end if 
   end if

   if (allocated(fix_weighted_for_element)) then
      if (allocated(fix_weighted_for_element_ini)) deallocate(fix_weighted_for_element_ini)
      allocate(fix_weighted_for_element_ini(size(fix_weighted_for_element)))
   end if

   if (allocated(fix_weighted_for_element_3ch)) then
      if (allocated(fix_weighted_for_element_3ch_ini)) deallocate(fix_weighted_for_element_3ch_ini)
      allocate(fix_weighted_for_element_3ch_ini(size(fix_weighted_for_element_3ch)))
   end if !weighted

   ! invisible chemical element set up ... I do it what
   int_local = nitems2(chemical_elements_invisible)
   linvisible = .false.
   if (int_local > 0) linvisible = .true.
   if (linvisible) then
      fix_no_of_elements_invisible = int_local

      if (allocated(fix_ch_elements_invisible)) deallocate (fix_ch_elements_invisible)
      allocate (fix_ch_elements_invisible(int_local))
      if (allocated(fix_Z_elements_invisible)) deallocate (fix_Z_elements_invisible)
      allocate (fix_Z_elements_invisible(int_local))
      if (allocated(fix_mass_elements_invisible)) deallocate (fix_mass_elements_invisible)
      allocate (fix_mass_elements_invisible(int_local))
      if (allocated(fix_type_to_periodic_invisible)) deallocate (fix_type_to_periodic_invisible)
      allocate (fix_type_to_periodic_invisible(int_local))


      icount = 0
      read (chemical_elements_invisible, *) fix_ch_elements_invisible(:)
      tmp = '-- list of invisible atoms --'//nwl// &
         ' id element   Z     mass'//nwl
      do int_local = 1, size(fix_ch_elements_invisible, dim=1)
         do i_p = 1, size_periodic_table
            if (periodic_table_element(i_p)%symbol == fix_ch_elements_invisible(int_local)) then
               icount = icount + 1
               fix_Z_elements_invisible(int_local) = periodic_table_element(i_p)%Z
               fix_mass_elements_invisible(int_local) = periodic_table_element(i_p)%mass
               fix_type_to_periodic_invisible(int_local) = i_p
            end if
         end do

         tmp = tmp//vtoa(int_local, 3)//vtoa(fix_type_to_periodic_invisible(int_local), 4)//fix_ch_elements_invisible(int_local)// &
            vtoa(fix_Z_elements_invisible(int_local))//vtoa(fix_mass_elements_invisible(int_local))//nwl

      end do
      call log_info(tmp)
   end if

   ! begin_species  defines type of interactions between species ...

   size_species_full = fix_no_of_elements**2
   size_species_half = fix_no_of_elements*(fix_no_of_elements+1)/2

   if (allocated(typ_species_half)) deallocate(typ_species_half)
   allocate(typ_species_half(size_species_half))

   if (allocated(typ_species_full)) deallocate(typ_species_full)
   allocate(typ_species_full(size_species_full))

   if (allocated(map_species_half)) deallocate(map_species_half)
   allocate(map_species_half(size_species_half,size_species_half))

   if (allocated(map_species_full)) deallocate(map_species_full)
   allocate(map_species_full(size_species_full,size_species_full))


   icount = 0
   do ii = 1, fix_no_of_elements
      icount = icount + 1
      map_species_half(ii, ii) = icount
      typ_species_half(icount)%type1 = ii
      typ_species_half(icount)%type2 = ii
      typ_species_half(icount)%Z1 = fix_Z_elements(ii)
      typ_species_half(icount)%Z2 = fix_Z_elements(ii)
      typ_species_half(icount)%ch1 = fix_ch_elements(ii)
      typ_species_half(icount)%ch2 = fix_ch_elements(ii)
   end do
   do ii = 1, fix_no_of_elements
      do jj = ii + 1, fix_no_of_elements
         icount = icount +1
         map_species_half(ii,jj) = icount
         map_species_half(jj,ii) = icount
         typ_species_half(icount)%type1 = ii
         typ_species_half(icount)%type2 = jj
         typ_species_half(icount)%Z1 = fix_Z_elements(ii)
         typ_species_half(icount)%Z2 = fix_Z_elements(jj)
         typ_species_half(icount)%ch1 = fix_ch_elements(ii)
         typ_species_half(icount)%ch2 = fix_ch_elements(jj)

      end do
   end do

   icount = 0
   do ii = 1, fix_no_of_elements
      do jj = 1, fix_no_of_elements
         icount = icount  +1
         map_species_full(ii, jj)  = icount
         typ_species_full(icount)%type1 = ii
         typ_species_full(icount)%type2 = jj
         typ_species_full(icount)%Z1 = fix_Z_elements(ii)
         typ_species_full(icount)%Z2 = fix_Z_elements(jj)
         typ_species_full(icount)%ch1 = fix_ch_elements(ii)
         typ_species_full(icount)%ch2 = fix_ch_elements(jj)
      end do
   end do
   ! end_species .......................

end subroutine fix_type_of_atoms




subroutine write_all_input

   use module_afs, only: n_rbf_afs, n_rbf, n_cheb, afs_type

   use ml_in_ndm_module, ONLY: ml_type,  rangml, &
      descriptor_type,  &
      nd_fingerprint,  &
      iread_ml, isave_ml, &
      toy_model, debug,  &
      n_g2_eta, n_g2_rs, n_g3_eta, n_g3_zeta, n_g3_lambda, &
      l_max, j_max, &
      alpha_soap, n_soap, lsoap, lsoap_fcut_wes, lsoap_diag, &
      weighted_auto,  &
      seed, &
      write_desc, write_desc_dump, read_desc_dump, &
      weighted, weighted_3ch, &
      marginal_likelihood, strict_behler,  &
      sign_stress, sign_stress_big_box, mld_fit_type, online_fit_batch_size, snap_class_constraints, &
      train_only,  r_cut_width_soap, &
      mtp_poly_min, mtp_poly_max,  &
      desc_forces, &
      eta_max_g2, eta_min_g2, rs_min_g2, rs_max_g2,  &
      lsoap_diag, lsoap_norm, lsoap_lnorm, nspecies_soap, &
      lambda_krr, min_lambda_krr, max_lambda_krr, n_values_lambda_krr,  &
      chemical_elements_invisible, mld_order, &
      polyc_n_poly, polyc_n_hermite,  &
      n_pca, classes_for_mcd,  rmat_dim,  dim_fourier_nbody, &
      length_fourier_nbody, delta_fourier_nbody, mld_regularization_type,  &
   !this is double for moment will remains only mld_
      mld_type_quadratic, &
      power_line, power_coeff_renorm, &
      write_design_matrix, write_test_design_matrix, &
      svd_rcond, type_of_loss, lmask, mask_file, &
      desc_file_format

   use module_bispectrum_so4, only: inv_r0_input, lbso4_diag
   use module_body_desc, only: bond_dist_transform, bond_beta, l_body_order, body_D_max, bond_dist_ann
   use module_kernel, only: kernel_type,  &
      length_kernel, sigma_kernel, write_kernel_matrix, &
      kernel_dump, np_kernel_ref, np_kernel_full, np_omega, &
      kernel_power, power_mcd, krff_type
   use module_nlinear, only: order_nlinear
   use k_cross_validation, only: kcross, n_kcross
   use def_kernels, ONLY: length_kse, sigma_kse
   use temporary_data_cov, only: dim_data
   use module_lbfgs_input, only: lbfgs_xtol, lbfgs_eps, &
      lbfgs_m_hess, lbfgs_max_steps, &
      lbfgs_print, lbfgs_gtol
   use module_db_poscar, only: iread_energy, ref_energy_per_element
   use module_so3, only: n_rbf_so3, radial_pow_so3
   use module_neigh_local, only: r_cut, r_cut_width, r_cut_in, r_cut_width_in,  type_fcut
   use module_ml_scalapack, only: scalapack_driver, debug_scalapack, nbr_predefined, nbc_predefined
   use module_zbl, only: zbl_potential, zbl_type, zbl_mode_default_k2b, zbl_mode_alone, r1_zbl, r2_zbl
   use module_variance, only : energy_variance, force_variance, stress_variance
   use module_Sigma_matrix, only: method_sigma
   use module_cur, only : cur_eps, cur_kval, cur_cval, cur_rval
   use module_optimization, only: optimize_weights_db, optimize_weights_chem, optimize_weights_L1, optimize_weights_L2, optimize_weights_Le, &
      optimize_ga_population, class_no_optimize_weights,  max_iter_optimize_weights, &
      factor_force_error, factor_energy_error, factor_stress_error
   use module_kernel_2b, only : activate_k2b, sigma_2b, delta_2b, np_radial_2b, r_cut_2b, r_cut_width_2b
   use module_chemical_species, only: fix_no_of_elements, chemical_elements,weight_per_element, weight_per_element_3ch, &
      img_num_ch, img_weighted, ws1, ws2, ws3, ws4
   use module_ftnbody, only : find_best_length, &
                             r_cut_ft2b, r_cut_width_ft2b, r_cut_ft3b, r_cut_width_ft3b, &
                             r_cut_ft4b, r_cut_width_ft4b, r_cut_ft5b, r_cut_width_ft5b, &
                             ftnbody_chem_mode, ftnbody_chem_rank, ftnbody_hash_channels, &
                             ftnbody_model, spip_degree_nbody
   use module_db_setup, only: selection_type, db_file, db_path, drop_short_dist
   use module_kernel_zetabody, only: dim_grid_zetabody,  dim_length_zetabody, dim_delta_zetabody, zetabody_order                          
   use module_ace_desc, only: l_ace_order, ace_lmax_list, ace_nmax_list, ace_lambda_list, ace_kmax_list, ace_numax, &
                              ace_radial_poly, ace_radial_chem, ace_gencg, & 
                              zetaace_order, dim_delta_zetaace, &
                              l_ace_set_rcut,  &
                              ace_rcut_in_list, ace_rcut_out_list, ace_rcut_width_in_list, ace_rcut_width_out_list, ace_chem, &
                              ace_npoints_spline, &
                              ace_chem_low_rank, ace_chem_low_rank_q, ace_chem_low_rank_niter, ace_chem_low_rank_lambda, &
                              ace_svd_randomized, ace_svd_randomized_oversample, ace_svd_randomized_power_iter
   use module_input_als_fit, only: als_nsteps, als_tol, als_ridge_k, als_rho, &
      als_alpha_method, als_precond_type, als_nnls_alpha, als_nnls_mode, als_nu_max, als_block_partition
   use module_tbind, only: tb_ham_ss, tb_ham_pp, tb_ham_sp, tb_ham_dd, &
      tb_nn_max, tb_kmax_list, tb_nmax_list, tb_lmax_list, tb_lambda_list, &
      tb_kmax_grid, tb_lmax_grid, &
      tb_rcut_in_list, tb_rcut_out_list, tb_rcut_width_in_list, tb_rcut_width_out_list, &
      tb_rcut_in_ij_list, tb_rcut_out_ij_list, tb_rcut_width_in_ij_list, tb_rcut_width_out_ij_list, &
      tb_model_lambda, tb_model_trace, tb_model_equivb, tb_power_trace_list, &
      tb_filter_type, tb_filter_order, tb_filter_ncenter, tb_filter_width, tb_filter_width_ratio, &
      tb_filter_use_first_moment, tb_filter_use_real_part, tb_filter_use_imag_part, &
      tb_filter_degree, tb_filter_spectral_padding, tb_filter_lmin_list, tb_filter_lmax_list, &
      tb_nrf, tb_nrg, tb_g_type, tb_svd, tbind_dim, TB_N_CHANNELS

   implicit none
   integer :: unumber, ii 


   if (rangml == 0) then
      open (newunit=unumber, file="old_input.ml", status='unknown')

      write (unumber, '("& input_ml ")')
      if (debug) then
         write (unumber, '("debug = .true. ")')
      else
         write (unumber, '("debug = .false. ")')
      end if
      write (unumber, '("ml_type = ",i5)') ml_type
      if (train_only) then
         write (unumber, '("train_only = .true. ")')
      else
         write (unumber, '("train_only = .false. ")')
      end if

      if (write_design_matrix) then
         write (unumber, '("write_design_matrix = .true. ")')
      else
         write (unumber, '("write_design_matrix = .false. ")')
      end if

      if (write_test_design_matrix) then
         write (unumber, '("write_test_design_matrix = .true. ")')
      else
         write (unumber, '("write_test_design_matrix = .false. ")')
      end if


      write (unumber, '("desc_file_format = ", i5)') desc_file_format

      write (unumber, '("!-------------------Scalapack Options-----------------------!")')

      if (scalapack_driver) then
         write (unumber, '("scalapack_driver = .true. ")')
      else
         write (unumber, '("scalapack_driver = .false. ")')
      end if
      if (debug_scalapack) then
         write (unumber, '("debug_scalapack = .true. ")')
      else
         write (unumber, '("debug_scalapack = .false. ")')
      end if
      write (unumber, '("npr_predefined = ", i5)') nbr_predefined
      write (unumber, '("npc_predefined = ", i5)') nbc_predefined
      write (unumber, '("!-------------------Function Fit form-----------------------!")')
      write (unumber, '("mld_order = ", i5)') mld_order
      write (unumber, '("mld_type_quadratic = ", i5)') mld_type_quadratic
      write (unumber, '("order_nlinear = ", i5)') order_nlinear
      write (unumber, '("polyc_n_poly = ", i5)') polyc_n_poly
      write (unumber, '("polyc_n_hermite = ", i5)') polyc_n_hermite

      write (unumber, '("!-----------------Algo Fit / regularization form--------------!")')
      write (unumber, '("mld_fit_type = ", i5)') mld_fit_type
      write (unumber, '("online_fit_batch_size = ", i5)') online_fit_batch_size
      write (unumber, '("mld_regularization_type = ", i5)') mld_regularization_type
      write (unumber, '("snap_class_constraints = """, a,"""")') trim(snap_class_constraints)
      write (unumber, '("lambda_krr = ", d22.15)') lambda_krr
      write (unumber, '("min_lambda_krr = ", d22.15)') min_lambda_krr
      write (unumber, '("max_lambda_krr = ", d22.15)') max_lambda_krr
      write (unumber, '("svd_rcond = ", d22.15)') svd_rcond
      write (unumber, '("n_values_lambda_krr = ", i6)') n_values_lambda_krr
      write(unumber,  '("type_of_loss = ", i6)') type_of_loss


      write (unumber, '("!-----------------ALS-Ridge (mld_fit_type=5)------------------!")')
      write (unumber, '("als_nsteps = ", i6)') als_nsteps
      write (unumber, '("als_tol = ", d22.15)') als_tol
      write (unumber, '("als_ridge_k = ", d22.15)') als_ridge_k
      write (unumber, '("als_rho = ", d22.15)') als_rho
      write (unumber, '("als_alpha_method = ", i6)') als_alpha_method
      write (unumber, '("als_precond_type = ", i6)') als_precond_type
      if (als_nnls_alpha) then
         write (unumber, '("als_nnls_alpha = .true. ")')
      else
         write (unumber, '("als_nnls_alpha = .false. ")')
      end if
      write (unumber, '("als_nnls_mode = ", i6)') als_nnls_mode
      write (unumber, '("als_nu_max = ", i6)') als_nu_max
      if (allocated(als_block_partition)) then
         do ii = 1, size(als_block_partition)
            write (unumber, '("als_block_partition(", i3, ") = ", i9)') ii, als_block_partition(ii)
         end do
      end if


      write (unumber, '("!---------------------------Kernel / -------------------------!")')

      ! some kernel stuffs.
      if (write_kernel_matrix) then
         write (unumber, '("write_kernel_matrix = .true. ")')
      else
         write (unumber, '("write_kernel_matrix = .false. ")')
      end if
      write (unumber, '("kernel_type = ", i6)') kernel_type
      write (unumber, '("kernel_rff = ", i6)') krff_type
      write (unumber, '("np_kernel_ref = ", i6)') np_kernel_ref
      write (unumber, '("np_kernel_full = ", i6)') np_kernel_full
      write (unumber, '("np_omega = ", i6)') np_omega
      write (unumber, '("kernel_dump = ", i6)') kernel_dump
      write (unumber, '("length_kse = ", d22.15)') length_kse
      write (unumber, '("sigma_kse = ", d22.15)') sigma_kse
      write (unumber, '("length_kernel = ", d22.15)') length_kernel
      write (unumber, '("sigma_kernel = ", d22.15)') sigma_kernel
      write (unumber, '("kernel_power = ", d22.15)') kernel_power
      write (unumber, '("power_mcd = ", d22.15)') power_mcd
      write (unumber, '("method_sigma = ", i6)') method_sigma
      write (unumber, '("cur_kval = ", i6)') cur_kval
      write (unumber, '("cur_cval = ", i6)') cur_cval
      write (unumber, '("cur_rval = ", i6)') cur_rval
      write (unumber, '("cur_eps = ", d22.15)') cur_eps

      if (activate_k2b) then
         write (unumber, '("activate_k2b = .true. ")')
      else
         write (unumber, '("activate_k2b = .false. ")')
      end if

      write (unumber, '("sigma_2b = ", d22.15)') sigma_2b
      write (unumber, '("delta_2b = ", d22.15)') delta_2b
      write (unumber, '("np_radial_2b = ", i6)') np_radial_2b
      write (unumber, '("r_cut_2b = ", d22.15)') r_cut_2b
      write (unumber, '("r_cut_width_2b = ", d22.15)') r_cut_width_2b

      write (unumber, '("!----------------------DB and Elements ------------------------!")')

      write (unumber, '("db_file = """, a, """")') trim(db_file)
      write (unumber, '("db_path = """, a, """")') trim(db_path)
      write (unumber, '("drop_short_dist = ", d22.15)') drop_short_dist
      write (unumber, '("selection_type = ", i5)') selection_type
      write (unumber, '("seed = ", i9)') seed
      if (weighted) then
         write (unumber, '("weighted = .true. ")')
      else
         write (unumber, '("weighted = .false. ")')
      end if

      if (weighted_auto) then
         write (unumber, '("weighted_auto = .true. ")')
      else
         write (unumber, '("weighted_auto = .false. ")')
      end if

      if (weighted_3ch) then
         write (unumber, '("weighted_3ch = .true. ")')
      else
         write (unumber, '("weighted_3ch = .false. ")')
      end if


      write (unumber, '("fix_no_of_elements = ", i5)') fix_no_of_elements
      write (unumber, '("iread_energy = ", i5)') iread_energy
      write (unumber, '("chemical_elements = """, a, """")') trim(chemical_elements)
      write (unumber, '("chemical_elements_invisible = """, a, """")') trim(chemical_elements_invisible)
      write (unumber, '("weight_per_element = """, a, """")') trim(weight_per_element)
      write (unumber, '("weight_per_element_3ch = """, a, """")') trim(weight_per_element_3ch)
      write (unumber, '("ref_energy_per_element = """, a, """")') trim(ref_energy_per_element)

      write (unumber, '("n_pca = ", i5)') n_pca
      write (unumber, '("classes_for_mcd = """, a, """")') trim(classes_for_mcd)
      write (unumber, '("sign_stress = ", d22.15)') sign_stress
      write (unumber, '("sign_stress_big_box = ", d22.15)') sign_stress_big_box
      if (zbl_potential) then
         write (unumber, '("zbl_potential = .true. ")')
      else
         write (unumber, '("zbl_potential = .false. ")')
      end if
      write (unumber, '("zbl_type = ", i6)') zbl_type
      write (unumber, '("r1_zbl = ", d22.15)') r1_zbl
      write (unumber, '("r2_zbl = ", d22.15)') r2_zbl


      if (force_variance) then
         write (unumber, '("force_variance = .true. ")')
      else
         write (unumber, '("force_variance = .false. ")')
      end if
      if (energy_variance) then
         write (unumber, '("energy_variance = .true. ")')
      else
         write (unumber, '("energy_variance = .false. ")')
      end if
      if (stress_variance) then
         write (unumber, '("stress_variance = .true. ")')
      else
         write (unumber, '("stress_variance = .false. ")')
      end if




      write (unumber, '("!------------------------Descriptors------------------------!")')

      write (unumber, '("descriptor_type = ", i6)') descriptor_type
      if (write_desc) then
         write (unumber, '("write_desc = .true. ")')
      else
         write (unumber, '("write_desc = .false. ")')
      end if
      if (lmask) then
         write (unumber, '("lmask = .true. ")')
      else
         write (unumber, '("lmask = .false. ")')
      end if

      if (mask_file) then
         write (unumber, '("mask_file = .true. ")')
      else
         write (unumber, '("mask_file = .false. ")')
      end if


      if (write_desc_dump) then
         write (unumber, '("write_desc_dump = .true. ")')
      else
         write (unumber, '("write_desc_dump = .false. ")')
      end if


      if (read_desc_dump) then
         write (unumber, '("read_desc_dump = .true. ")')
      else
         write (unumber, '("read_desc_dump = .false. ")')
      end if

      if (desc_forces) then
         write (unumber, '("desc_forces = .true. ")')
      else
         write (unumber, '("desc_forces = .false. ")')
      end if

      write (unumber, '("r_cut = ", d22.15)') r_cut
      write (unumber, '("r_cut_in = ", d22.15)') r_cut_in
      write (unumber, '("r_cut_width_in = ", d22.15)') r_cut_width_in
      write (unumber, '("r_cut_width = ", d22.15)') r_cut_width
      write (unumber, '("type_fcut = ", i6)') type_fcut

      write (unumber, '("!---G2----")')

      write (unumber, '("n_g2_eta = ", i6)') n_g2_eta
      write (unumber, '("n_g2_rs = ", i6)') n_g2_rs
      write (unumber, '("eta_max_g2 = ", d22.15)') eta_max_g2
      write (unumber, '("eta_min_g2 = ", d22.15)') eta_min_g2
      write (unumber, '("rs_max_g2 = ", d22.15)') rs_max_g2
      write (unumber, '("rs_min_g2 = ", d22.15)') rs_min_g2

      write (unumber, '("!---G3----")')

      write (unumber, '("n_g3_eta = ", i6)') n_g3_eta
      write (unumber, '("n_g3_zeta = ", i6)') n_g3_zeta
      write (unumber, '("n_g3_lambda = ", i6)') n_g3_lambda

      write (unumber, '("!---BehPar----")')
      if (strict_behler) then
         write (unumber, '("strict_behler = .true. ")')
      else
         write (unumber, '("strict_behler = .false. ")')
      end if


      write (unumber, '("!---AFS----")')
      write (unumber, '("afs_type = ", i6)') afs_type
      write (unumber, '("n_rbf = ", i6)') n_rbf
      write (unumber, '("n_rbf_afs = ", i6)') n_rbf_afs
      write (unumber, '("n_rbf_so3 = ", i6)') n_rbf_so3
      write (unumber, '("n_cheb = ", i6)') n_cheb


      write (unumber, '("!---MTP----")')
      write (unumber, '("mtp_poly_min = ", i6)') mtp_poly_min
      write (unumber, '("mtp_poly_max = ", i6)') mtp_poly_max


      write (unumber, '("!---SO3----")')
      write (unumber, '("l_max = ", i6)') l_max
      write (unumber, '("radial_pow_so4 = ", i6)') radial_pow_so3


      write (unumber, '("!---SO4/BS04----")')
      write (unumber, '("j_max = ", f15.10)') j_max
      write (unumber, '("inv_r0_input = ", d22.15)') inv_r0_input
      if (lbso4_diag) then
         write (unumber, '("lbso4_diag = .true. ")')
      else
         write (unumber, '("lbso4_diag = .false. ")')
      end if



      write (unumber, '("!---SOAP----")')
      write (unumber, '("alpha_soap = ", f15.10)') alpha_soap
      write (unumber, '("n_soap = ", i6)') n_soap
      write (unumber, '("r_cut_width_soap = ", d22.15)') r_cut_width_soap
      write (unumber, '("r_cut_width      = ", d22.15)') r_cut_width
      if (lsoap_diag) then
         write (unumber, '("lsoap_diag = .true. ")')
      else
         write (unumber, '("lsoap_diag = .false. ")')
      end if

      if (lsoap_norm) then
         write (unumber, '("lsoap_norm = .true. ")')
      else
         write (unumber, '("lsoap_norm = .false. ")')
      end if

      if (lsoap_lnorm) then
         write (unumber, '("lsoap_lnorm = .true. ")')
      else
         write (unumber, '("lsoap_lnorm = .false. ")')
      end if
      write (unumber, '("nspecies_soap = ", i6)') nspecies_soap

      if (lsoap) then
         write (unumber, '("lsoap = .true. ")')
      else
         write (unumber, '("lsoap = .false. ")')
      end if

      if (lsoap_fcut_wes) then
         write (unumber, '("lsoap_fcut_wes = .true. ")')
      else
         write (unumber, '("lsoap_fcut_wes = .false. ")')
      end if

      write (unumber, '("!---BODY----")')
      write (unumber, '("body_D_max(1) = ", i6)') body_D_max(1)
      write (unumber, '("body_D_max(2) = ", i6)') body_D_max(2)
      write (unumber, '("body_D_max(3) = ", i6)') body_D_max(3)
      write (unumber, '("body_D_max(4) = ", i6)') body_D_max(4)
      write (unumber, '("body_D_max(5) = ", i6)') body_D_max(5)

      if (l_body_order(1)) then
         write (unumber, '("l_body_order(1) = .true. ")')
      else
         write (unumber, '("l_body_order(1) = .false. ")')
      end if

      if (l_body_order(2)) then
         write (unumber, '("l_body_order(2) = .true. ")')
         write (unumber, '("dim_length_zetabody(2) = """, a, """")') trim(dim_length_zetabody(2))
         write (unumber, '("dim_grid_zetabody(2) = """, a, """")') trim(dim_grid_zetabody(2))
         write (unumber, '("dim_delta_zetabody(2) = """, a, """")') trim(dim_delta_zetabody(2))

      else
         write (unumber, '("l_body_order(2) = .false. ")')
      end if


      if (l_body_order(3)) then
         write (unumber, '("l_body_order(3) = .true. ")')
         write (unumber, '("dim_length_zetabody(3) = """, a, """")') trim(dim_length_zetabody(3))
         write (unumber, '("dim_grid_zetabody(3) = """, a, """")') trim(dim_grid_zetabody(3))
         write (unumber, '("dim_delta_zetabody(3) = """, a, """")') trim(dim_delta_zetabody(3))
      else
         write (unumber, '("l_body_order(3) = .false. ")')
      end if


      if (l_body_order(4)) then
         write (unumber, '("l_body_order(4) = .true. ")')
      else
         write (unumber, '("l_body_order(4) = .false. ")')
      end if


      if (l_body_order(5)) then
         write (unumber, '("l_body_order(5) = .true. ")')
      else
         write (unumber, '("l_body_order(5) = .false. ")')
      end if

      write (unumber, '("bond_dist_transform = ", i6)') bond_dist_transform

      write (unumber, '("bond_beta = ", d22.15)') bond_beta
      write (unumber, '("bond_dist_ann = ", d22.15)') bond_dist_ann

      write (unumber, '("zetabody_order = ", i6)') zetabody_order


      write (unumber, '("!----ACE----")')
      write (unumber, '("ace_chem = ", i6)') ace_chem
      write (unumber, '("ace_numax = ", i6)') ace_numax
      ! write (unumber, '("acekmax = ", i6)') acekmax
      write (unumber, '("ace_radial_poly = ", i6)') ace_radial_poly
      write (unumber, '("ace_radial_chem = ", i6)') ace_radial_chem
      write (unumber, '("ace_gencg = ", i6)') ace_gencg
      write (unumber, '("ace_npoints_spline = ", i6)') ace_npoints_spline
      write (unumber, '("ace_chem_low_rank = ", i6)') ace_chem_low_rank
      write (unumber, '("ace_chem_low_rank_q = ", i6)') ace_chem_low_rank_q
      write (unumber, '("ace_chem_low_rank_niter = ", i6)') ace_chem_low_rank_niter
      write (unumber, '("ace_chem_low_rank_lambda = ", ES12.4)') ace_chem_low_rank_lambda
      write (unumber, '("ace_svd_randomized = ", i6)') ace_svd_randomized
      write (unumber, '("ace_svd_randomized_oversample = ", i6)') ace_svd_randomized_oversample
      write (unumber, '("ace_svd_randomized_power_iter = ", i6)') ace_svd_randomized_power_iter
      do ii = 1, ace_numax
         if (l_ace_order(ii)) write (unumber, '("l_ace_order(", i6, ") = .true. ")') ii
      end do
      write (unumber, '("ace_nmax_list = """, a,"""")') ace_nmax_list
      write (unumber, '("ace_kmax_list = """, a,"""")') ace_kmax_list
      write (unumber, '("ace_lmax_list = """, a,"""")') ace_lmax_list
      write (unumber, '("ace_lambda_list =""", a,"""")') ace_lambda_list

      write (unumber, '("zetaace_order = ", i6)') zetaace_order

      if (zetaace_order > 1) then 
        do ii = 1, zetaace_order
          if (l_ace_order(ii)) then
            write(unumber, '("dim_delta_zetaace(", i6, ") = """, a, """")') ii, trim(dim_delta_zetaace(ii))
          end if 
        end do 
      end if 
      

      !write (unumber, '("ace_lambda_list(", i6, ") =""", a,"""")') ii, ace_lambda_list(ii)

      if (l_ace_set_rcut) then
         write (unumber, '("l_ace_set_rcut = .true. ")')
      else
         write (unumber, '("l_ace_set_rcut = .false. ")')
      end if


      write (unumber, '("ace_rcut_in_list =""", a,"""")') ace_rcut_in_list
      write (unumber, '("ace_rcut_width_in_list =""", a,"""")') ace_rcut_width_in_list
      write (unumber, '("ace_rcut_out_list =""", a,"""")') ace_rcut_out_list
      write (unumber, '("ace_rcut_width_out_list =""", a,"""")') ace_rcut_width_out_list
      write (unumber, '("ace_lambda_list =""", a,"""")') ace_lambda_list 
      







      write (unumber, '("!----MiLaDy----")')
      write (unumber, '("rmat_dim = ", i6)') rmat_dim
      write (unumber, '("img_num_ch = ", i6)') img_num_ch
      write (unumber, '("power_line = ", d22.15)') power_line
      write (unumber, '("power_coeff_renorm = ", d22.15)') power_coeff_renorm
      if (img_weighted) then
         write (unumber, '("img_weighted = .true. ")')
      else
         write (unumber, '("img_weighted = .false. ")')
      end if
      write (unumber, '("ws1 = """, a, """")') trim(ws1)
      write (unumber, '("ws2 = """, a, """")') trim(ws2)
      write (unumber, '("ws3 = """, a, """")') trim(ws3)
      write (unumber, '("ws4 = """, a, """")') trim(ws4)



      write (unumber, '("!---DIRECT NBODY----")')
      do ii = 1, size(l_body_order)
         if (l_body_order(ii)) write (unumber, '(" body  =", i6, "dim_grid_zetabody = """, a, """")') ii, trim(dim_grid_zetabody(ii))
      end do 
      ! if (descriptor_type == descriptor_ftnbody) then
      write (unumber, '("!---FOURIER NBODY----")')
      write (unumber, '("dim_fourier_nbody = """, a, """")') trim(dim_fourier_nbody)
      write (unumber, '("length_fourier_nbody = """, a, """")') trim(length_fourier_nbody)
      write (unumber, '("delta_fourier_nbody = """, a, """")') trim(delta_fourier_nbody)
      write (unumber, '("ftnbody_chem_mode = ", i6)') ftnbody_chem_mode
      write (unumber, '("ftnbody_chem_rank = """, a, """")') trim(ftnbody_chem_rank)
      write (unumber, '("ftnbody_hash_channels = """, a, """")') trim(ftnbody_hash_channels)
      write (unumber, '("ftnbody_model = """, a, """")') trim(ftnbody_model)
      write (unumber, '("spip_degree_nbody = """, a, """")') trim(spip_degree_nbody)
    if (find_best_length) then
      write (unumber, '("find_best_length = .true. ")')
    else
      write (unumber, '("find_best_length = .false. ")')
    end if
    write (unumber, '("r_cut_ft2b = ", f15.10)') r_cut_ft2b
    write (unumber, '("r_cut_width_ft2b = ", f15.10)') r_cut_width_ft2b
    write (unumber, '("r_cut_ft3b = ", f15.10)') r_cut_ft3b
    write (unumber, '("r_cut_width_ft3b = ", f15.10)') r_cut_width_ft3b
    write (unumber, '("r_cut_ft4b = ", f15.10)') r_cut_ft4b
    write (unumber, '("r_cut_width_ft4b = ", f15.10)') r_cut_width_ft4b
    write (unumber, '("r_cut_ft5b = ", f15.10)') r_cut_ft5b
    write (unumber, '("r_cut_width_ft5b = ", f15.10)') r_cut_width_ft5b
      ! end if


      write (unumber, '("!----------------------Genetical algo-----------------------!")')


      if (optimize_weights_db) then
         write (unumber, '("optimize_weights_db = .true. ")')
      else
         write (unumber, '("optimize_weights_db = .false. ")')
      end if

      if (optimize_weights_chem) then
         write (unumber, '("optimize_weights_chem = .true. ")')
      else
         write (unumber, '("optimize_weights_chem = .false. ")')
      end if


      if (optimize_weights_L1) then
         write (unumber, '("optimize_weights_L1 = .true. ")')
      else
         write (unumber, '("optimize_weights_L1 = .false. ")')
      end if

      if (optimize_weights_L2) then
         write (unumber, '("optimize_weights_L2 = .true. ")')
      else
         write (unumber, '("optimize_weights_L2 = .false. ")')
      end if
      if (optimize_weights_Le) then
         write (unumber, '("optimize_weights_Le = .true. ")')
      else
         write (unumber, '("optimize_weights_Le = .false. ")')
      end if

      write (unumber, '("optimize_ga_population = ", i6)') optimize_ga_population
      write (unumber, '("class_no_optimize_weights = """, a, """")') trim(class_no_optimize_weights)
      write (unumber, '("max_iter_optimize_weights = ", i6)') max_iter_optimize_weights
      ! weights optimization using genetic algorithms ...

      write (unumber, '("factor_energy_error = ", d22.15)') factor_energy_error
      write (unumber, '("factor_force_error = ", d22.15)') factor_force_error
      write (unumber, '("factor_stress_error = ", d22.15)') factor_stress_error


      write (unumber, '("!-----------------------TB DESCRIPTOR-----------------------!")')
      write (unumber, '("descriptor_type = ", i6)') descriptor_type
      if (tb_ham_ss) then
         write (unumber, '("tb_ham_ss = .true. ")')
      else
         write (unumber, '("tb_ham_ss = .false. ")')
      end if
      if (tb_ham_pp) then
         write (unumber, '("tb_ham_pp = .true. ")')
      else
         write (unumber, '("tb_ham_pp = .false. ")')
      end if
      if (tb_ham_sp) then
         write (unumber, '("tb_ham_sp = .true. ")')
      else
         write (unumber, '("tb_ham_sp = .false. ")')
      end if
      if (tb_ham_dd) then
         write (unumber, '("tb_ham_dd = .true. ")')
      else
         write (unumber, '("tb_ham_dd = .false. ")')
      end if
      write (unumber, '("tb_nn_max = ", 4(i6, ","))') tb_nn_max(1:TB_N_CHANNELS)
      write (unumber, '("tb_kmax_list = """, a, """")') trim(tb_kmax_list)
      write (unumber, '("tb_nmax_list = """, a, """")') trim(tb_nmax_list)
      write (unumber, '("tb_lmax_list = """, a, """")') trim(tb_lmax_list)
      write (unumber, '("tb_lambda_list = """, a, """")') trim(tb_lambda_list)
      write (unumber, '("tb_kmax_grid = """, a, """")') trim(tb_kmax_grid)
      write (unumber, '("tb_lmax_grid = """, a, """")') trim(tb_lmax_grid)
      write (unumber, '("tb_rcut_in_list = """, a, """")') trim(tb_rcut_in_list)
      write (unumber, '("tb_rcut_out_list = """, a, """")') trim(tb_rcut_out_list)
      write (unumber, '("tb_rcut_width_in_list = """, a, """")') trim(tb_rcut_width_in_list)
      write (unumber, '("tb_rcut_width_out_list = """, a, """")') trim(tb_rcut_width_out_list)
      write (unumber, '("tb_rcut_in_ij_list = """, a, """")') trim(tb_rcut_in_ij_list)
      write (unumber, '("tb_rcut_out_ij_list = """, a, """")') trim(tb_rcut_out_ij_list)
      write (unumber, '("tb_rcut_width_in_ij_list = """, a, """")') trim(tb_rcut_width_in_ij_list)
      write (unumber, '("tb_rcut_width_out_ij_list = """, a, """")') trim(tb_rcut_width_out_ij_list)
      write (unumber, '("tb_g_type = ", i6)') tb_g_type
      write (unumber, '("tb_svd = ", i6)') tb_svd
      if (tb_model_lambda) then
         write (unumber, '("tb_model_lambda = .true. ")')
      else
         write (unumber, '("tb_model_lambda = .false. ")')
      end if
      if (tb_model_trace) then
         write (unumber, '("tb_model_trace = .true. ")')
      else
         write (unumber, '("tb_model_trace = .false. ")')
      end if
      ! tbind_model is intentionally not re-emitted: the resolved tb_model_*
      ! logicals above carry the full model selection on restart.
      if (tb_model_equivb) then
         write (unumber, '("tb_model_equivb = .true. ")')
      else
         write (unumber, '("tb_model_equivb = .false. ")')
      end if
      write (unumber, '("tb_power_trace_list = """, a, """")') trim(tb_power_trace_list)
      write (unumber, '("tb_filter_type = ''", a, "''")') trim(tb_filter_type)
      write (unumber, '("tb_filter_order = ", i6)') tb_filter_order
      write (unumber, '("tb_filter_ncenter = ", i6)') tb_filter_ncenter
      write (unumber, '("tb_filter_width = ", d22.15)') tb_filter_width
      write (unumber, '("tb_filter_width_ratio = ", d22.15)') tb_filter_width_ratio
      if (tb_filter_use_first_moment) then
         write (unumber, '("tb_filter_use_first_moment = .true. ")')
      else
         write (unumber, '("tb_filter_use_first_moment = .false. ")')
      end if
      if (tb_filter_use_real_part) then
         write (unumber, '("tb_filter_use_real_part = .true. ")')
      else
         write (unumber, '("tb_filter_use_real_part = .false. ")')
      end if
      if (tb_filter_use_imag_part) then
         write (unumber, '("tb_filter_use_imag_part = .true. ")')
      else
         write (unumber, '("tb_filter_use_imag_part = .false. ")')
      end if
      write (unumber, '("tb_filter_degree = ", i6)') tb_filter_degree
      write (unumber, '("tb_filter_spectral_padding = ", d22.15)') tb_filter_spectral_padding
      write (unumber, '("tb_filter_lmin_list = """, a, """")') trim(tb_filter_lmin_list)
      write (unumber, '("tb_filter_lmax_list = """, a, """")') trim(tb_filter_lmax_list)
      write (unumber, '("tb_nrf = ", i6)') tb_nrf
      write (unumber, '("tb_nrg = ", i6)') tb_nrg
      write (unumber, '("tbind_dim = ", i6)') tbind_dim


      write (unumber, '("!--------------------------LBFGS----------------------------!")')

      write (unumber, '("lbfgs_m_hess = ", i6)') lbfgs_m_hess
      write (unumber, '("lbfgs_max_steps = ", i6)') lbfgs_max_steps
      write (unumber, '("lbfgs_print(1) = ", i6)') lbfgs_print(1)
      write (unumber, '("lbfgs_print(2) = ", i6)') lbfgs_print(2)
      write (unumber, '("lbfgs_eps = ", d22.15)') lbfgs_eps
      write (unumber, '("lbfgs_gtol = ", d22.15)') lbfgs_gtol
      write (unumber, '("lbfgs_xtol = ", d22.15)') lbfgs_xtol


      write (unumber, '("!-----------------------OLD UNUSED----------------------------!")')

      write (unumber, '("iread_ml = ", i6)') iread_ml
      write (unumber, '("isave_ml = ", i6)') isave_ml
      if (toy_model) then
         write (unumber, '("toy_model = .true. ")')
      else
         write (unumber, '("toy_model = .false. ")')
      end if


      if (kcross) then
         write (unumber, '("kcross = .true. ")')
      else
         write (unumber, '("kcross = .false. ")')
      end if

      if (marginal_likelihood) then
         write (unumber, '("marginal_likelihood = .true. ")')
      else
         write (unumber, '("marginal_likelihood = .false. ")')
      end if

      write (unumber, '("n_kcross = ", i9)') n_kcross
      write (unumber, '("dim_data = ", i9)') dim_data
      write (unumber, '("nd_fingerprint = ", i9)') nd_fingerprint
      write (unumber, '("max_data = ", i6)') n_kcross


      write (unumber, '("n_kcross = ", i6)') n_kcross



      write (unumber, '("& end ")')
      close (unumber, status='keep')

   end if

end subroutine write_all_input


subroutine fix_double_input_for_integer (name_true, ival_true_read, ival_true_default, &
                             name_fake, ival_fake_read, ival_fake_default, idefault)
   use mld_logger
   use mld_mpi
   implicit none
   integer, intent(in) :: ival_true_default, ival_fake_read, ival_fake_default
   integer, intent(inout) :: ival_true_read
   integer, intent(in) :: idefault
   character (len=*) :: name_true, name_fake
   integer, parameter :: imagic = -777
   integer :: izozo 

   izozo = ival_true_default
   izozo = ival_fake_default
   ! if both are changed ...
   if ((ival_true_read /= imagic ).and.(ival_fake_read /= imagic)) then
      if (ival_true_read /= ival_fake_read) then
         call mld_mpi_abort("read_ml :" //name_true//" and "//name_fake//" are incompatible. Choose only one! Or put them equal "//vtoa(ival_true_read)//vtoa(ival_fake_read))
      end if
   end if

   if ((ival_fake_read /= imagic ).and. (ival_true_read /= imagic )) then
      if (ival_true_read /= ival_fake_read) then
         call mld_mpi_abort("read_ml :" //name_true//" and "//name_fake//" are incompatible.  Choose only one! Or put them equal "//vtoa(ival_true_read)//vtoa(ival_fake_read))
      end if
   end if

   if ((ival_fake_read == imagic ).and. (ival_true_read == imagic )) then
      ival_true_read = idefault
   end if


   if ((ival_true_read /= imagic ).and.(ival_fake_read == imagic)) then
      ival_true_read = ival_true_read
   end if

   if ((ival_true_read ==imagic ).and.(ival_fake_read /= imagic)) then
      ival_true_read = ival_fake_read
   end if
end subroutine fix_double_input_for_integer

!call fix_double_input_for_integer(" mld_type_quadratic ", &
!                                     mld_type_quadratic, mld_type_quadratic_copy, &
!                                    " snap_type_quadratic ", &
!                                    snap_type_quadratic, snap_type_quadratic_copy, &
!                                    mld_type_quadratic_default)

subroutine fix_double_input_for_logical (name_true, lval_true_read, lval_true_default, &
   name_fake, lval_fake_read, lval_fake_default, ldefault)
   use mld_logger
   use mld_mpi
   implicit none
   logical, intent(inout) :: lval_true_default, lval_fake_read, lval_fake_default
   logical, intent(inout) :: lval_true_read
   logical, intent(in) :: ldefault
   character (len=*) :: name_true, name_fake
   logical :: lzozo 

   
   lzozo = lval_fake_default
   lzozo = lval_true_default
! if both are changed ...
   if ((lval_true_read .neqv. ldefault ).and.(lval_fake_read .neqv. ldefault)) then
      lval_true_read = lval_true_read
      lval_fake_read = lval_true_read
      call log_warning("read_ml :" //name_true//" and "//name_fake//" are in the same input at different values")
      call log_warning("read_ml :" //name_fake//" was removed form input. Here was initialized with "//name_true//" value ")
   end if

   if ((lval_true_read .neqv. ldefault ).and.(lval_fake_read .neqv. ldefault)) then
      lval_true_read = lval_true_read
   end if

   if ((lval_true_read .eqv. ldefault ).and.(lval_fake_read .neqv. ldefault)) then
      lval_true_read = ldefault
      lval_fake_read = ldefault
      call log_warning("read_ml :" //name_fake//" is present in the input list. Please use " // name_true// " instead ")
      call log_warning("read_ml :" //name_fake//" was removed form input list. Here was initialized with "//name_true//" default value ")
   end if


end subroutine fix_double_input_for_logical
