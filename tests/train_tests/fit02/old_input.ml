& input_ml 
debug = .false. 
ml_type =     0
train_only = .false. 
write_design_matrix = .false. 
write_test_design_matrix = .false. 
desc_file_format =     1
!-------------------Scalapack Options-----------------------!
scalapack_driver = .false. 
debug_scalapack = .false. 
npr_predefined =    32
npc_predefined =    32
!-------------------Function Fit form-----------------------!
mld_order =     1
mld_type_quadratic =     1
order_nlinear =     2
polyc_n_poly =     2
polyc_n_hermite =     1
!-----------------Algo Fit / regularization form--------------!
mld_fit_type =     4
mld_regularization_type =     0
snap_class_constraints = "06"
lambda_krr =  0.100000000000000D-07
min_lambda_krr =  0.100000000000000D-08
max_lambda_krr =  0.100000000000000D+02
svd_rcond = -0.100000000000000D+01
n_values_lambda_krr =     20
type_of_loss =      1
!-----------------ALS-Ridge (mld_fit_type=5)------------------!
als_nsteps =     10
als_tol = -0.100000000000000D+01
als_ridge_k =  0.100000000000000D-02
als_rho =  0.100000000000000D-05
als_alpha_method =      1
als_precond_type =      0
als_nnls_alpha = .false. 
als_nnls_mode =      0
als_nu_max =      1
!---------------------------Kernel / -------------------------!
write_kernel_matrix = .false. 
kernel_type =      1
kernel_rff =      1
np_kernel_ref =    200
np_kernel_full =    800
np_omega =      4
kernel_dump =      1
length_kse =  0.100000000000000D+01
sigma_kse =  0.100100000000000D+03
length_kernel =  0.707106781186548D+00
sigma_kernel =  0.000000000000000D+00
kernel_power =  0.200000000000000D+01
power_mcd =  0.500000000000000D-01
method_sigma =      1
cur_kval =     -1
cur_cval =     -1
cur_rval =     -1
cur_eps =  0.100000000000000D+01
activate_k2b = .false. 
sigma_2b =  0.500000000000000D-01
delta_2b =  0.100000000000000D+01
np_radial_2b =     30
r_cut_2b = -0.100000000000000D+01
r_cut_width_2b =  0.500000000000000D+00
!----------------------DB and Elements ------------------------!
db_file = "db_model.in"
db_path = "DB/"
selection_type =     1
seed =        37
weighted = .false. 
weighted_auto = .false. 
weighted_3ch = .false. 
fix_no_of_elements =     1
iread_energy =     2
chemical_elements = " W"
chemical_elements_invisible = ""
weight_per_element = "1.d0"
weight_per_element_3ch = "1.d0"
ref_energy_per_element = "0.d0"
n_pca =     3
classes_for_mcd = "07 08 12"
sign_stress =  0.100000000000000D+01
sign_stress_big_box =  0.100000000000000D+01
zbl_potential = .false. 
zbl_type =      1
r1_zbl =  0.100000000000000D+01
r2_zbl =  0.220000000000000D+01
force_variance = .false. 
energy_variance = .false. 
stress_variance = .false. 
!------------------------Descriptors------------------------!
descriptor_type =      4
write_desc = .false. 
lmask = .false. 
mask_file = .false. 
write_desc_dump = .false. 
read_desc_dump = .false. 
desc_forces = .true. 
r_cut =  0.470000000000000D+01
r_cut_in = -0.120000000000000D+01
r_cut_width_in =  0.400000000000000D+00
r_cut_width =  0.500000000000000D+00
type_fcut =      2
!---G2----
n_g2_eta =     29
n_g2_rs =      1
eta_max_g2 =  0.200000000000000D+01
eta_min_g2 =  0.200000000000000D+00
rs_max_g2 =  0.000000000000000D+00
rs_min_g2 =  0.000000000000000D+00
!---G3----
n_g3_eta =      7
n_g3_zeta =      3
n_g3_lambda =      2
!---BehPar----
strict_behler = .true. 
!---AFS----
afs_type =      1
n_rbf =     11
n_rbf_afs =     11
n_rbf_so3 =     11
n_cheb =     11
!---MTP----
mtp_poly_min =      2
mtp_poly_max =      7
!---SO3----
l_max =     11
radial_pow_so4 =      1
!---SO4/BS04----
j_max =    2.0000000000
inv_r0_input =  0.993633802276324D+00
lbso4_diag = .false. 
!---SOAP----
alpha_soap =    3.0000000000
n_soap =     11
r_cut_width_soap =  0.500000000000000D+00
r_cut_width      =  0.500000000000000D+00
lsoap_diag = .false. 
lsoap_norm = .true. 
lsoap_lnorm = .true. 
nspecies_soap =      1
lsoap = .false. 
lsoap_fcut_wes = .false. 
!---BODY----
body_D_max(1) =      1
body_D_max(2) =      2
body_D_max(3) =      3
body_D_max(4) =      4
body_D_max(5) =      5
l_body_order(1) = .true. 
l_body_order(2) = .true. 
dim_length_zetabody(2) = " 0.5d0"
dim_grid_zetabody(2) = " 30"
dim_delta_zetabody(2) = " 0.1d0 1.d0 10.d0"
l_body_order(3) = .true. 
dim_length_zetabody(3) = " 0.5d0 0.5d0 0.5d0"
dim_grid_zetabody(3) = " 10 10 20"
dim_delta_zetabody(3) = " 0.1d0 1.d0 10.d0"
l_body_order(4) = .true. 
l_body_order(5) = .false. 
bond_dist_transform =      3
bond_beta =  0.200000000000000D+01
bond_dist_ann =  0.100000000000000D+01
zetabody_order =      1
!----ACE----
ace_chem =      1
ace_numax =      3
ace_radial_poly =      2
ace_radial_chem =      1
ace_gencg =      2
ace_npoints_spline =   1000
ace_chem_low_rank =      0
ace_chem_low_rank_q =      4
ace_chem_low_rank_niter =     50
ace_chem_low_rank_lambda =   1.0000E-06
ace_nmax_list = " 5 3 1                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  "
ace_kmax_list = " 3 3 3                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  "
ace_lmax_list = " 0 4 2                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  "
ace_lambda_list =" 2.d0 2.d0 2.d0                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         "
zetaace_order =      1
l_ace_set_rcut = .false. 
ace_rcut_in_list =" 1.2d0 1.2d0 1.2d0                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      "
ace_rcut_width_in_list =" 0.4d0 0.4d0 0.4d0                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      "
ace_rcut_out_list =" 5.d0 5.d0 5.d0                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         "
ace_rcut_width_out_list =" 0.5d0 0.5d0 0.5d0                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      "
ace_lambda_list =" 2.d0 2.d0 2.d0                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         "
!----MiLaDy----
rmat_dim =     50
img_num_ch =      1
power_line =  0.200000000000000D+01
power_coeff_renorm =  0.333000000000000D+00
img_weighted = .false. 
ws1 = "1.d0"
ws2 = "1.d0"
ws3 = "1.d0"
ws4 = "1.d0"
!---DIRECT NBODY----
 body  =     1dim_grid_zetabody = "                                                                                                    "
 body  =     2dim_grid_zetabody = " 30"
 body  =     3dim_grid_zetabody = " 10 10 20"
 body  =     4dim_grid_zetabody = "                                                                                                    "
!---FOURIER NBODY----
dim_fourier_nbody = " 0 0 2000 0 0"
length_fourier_nbody = " 0.7d0 0.7d0 0.7d0 0.7d0 0.7d0"
delta_fourier_nbody = " 1.d0 1.d0 1.d0 1.d0 1.d0"
find_best_length = .true. 
r_cut_ft2b =    5.0000000000
r_cut_width_ft2b =    0.5000000000
r_cut_ft3b =    5.0000000000
r_cut_width_ft3b =    0.5000000000
r_cut_ft4b =    5.0000000000
r_cut_width_ft4b =    0.5000000000
r_cut_ft5b =    5.0000000000
r_cut_width_ft5b =    0.5000000000
!----------------------Genetical algo-----------------------!
optimize_weights_db = .false. 
optimize_weights_chem = .false. 
optimize_weights_L1 = .false. 
optimize_weights_L2 = .false. 
optimize_weights_Le = .true. 
optimize_ga_population =     80
class_no_optimize_weights = "02"
max_iter_optimize_weights =     20
factor_energy_error =  0.100000000000000D+01
factor_force_error =  0.100000000000000D+01
factor_stress_error =  0.100000000000000D+01
!-----------------------TB DESCRIPTOR-----------------------!
descriptor_type =      4
tb_ham_ss = .false. 
tb_ham_pp = .false. 
tb_ham_sp = .false. 
tb_ham_dd = .false. 
tb_nn_max =     -1,    -1,    -1,    -1,
tb_kmax_list = "8 8 8 8"
tb_nmax_list = "14 4 2 1"
tb_lmax_list = "3 2 2 1"
tb_lambda_list = "3.0 3.0 3.0 3.0"
tb_kmax_grid = "-1 -1 -1 -1"
tb_lmax_grid = "-1 -1 -1 -1"
tb_rcut_in_list = "1.2d0 1.2d0 1.2d0 1.2d0"
tb_rcut_out_list = "5.d0 5.d0 5.d0 5.d0"
tb_rcut_width_in_list = "0.4d0 0.4d0 0.4d0 0.4d0"
tb_rcut_width_out_list = "0.5d0 0.5d0 0.5d0 0.5d0"
tb_g_type =      1
tb_svd =      2
tb_model_lambda = .true. 
tb_model_trace = .false. 
tb_power_trace_list = "4 4 4 4"
tb_nrf =      1
tb_nrg =      1
tbind_dim =      0
!--------------------------LBFGS----------------------------!
lbfgs_m_hess =     40
lbfgs_max_steps =  20000
lbfgs_print(1) =    100
lbfgs_print(2) =      0
lbfgs_eps =  0.100000000000000D-01
lbfgs_gtol =  0.400000000000000D+00
lbfgs_xtol =  0.100000000000000D-15
!-----------------------OLD UNUSED----------------------------!
iread_ml =      0
isave_ml =      0
toy_model = .false. 
kcross = .false. 
marginal_likelihood = .false. 
n_kcross =         3
dim_data =      4000
nd_fingerprint =         3
max_data =      3
n_kcross =      3
& end 
