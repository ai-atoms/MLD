& input_ml 
debug = .false. 
ml_type =    -1
train_only = .true. 
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
polyc_n_hermite =     2
!-----------------Algo Fit / regularization form--------------!
mld_fit_type =     4
mld_regularization_type =     0
snap_class_constraints = "02"
lambda_krr = -0.100000000000000D+01
min_lambda_krr =  0.100000000000000D-09
max_lambda_krr =  0.100000000000000D+11
svd_rcond = -0.100000000000000D+01
n_values_lambda_krr =     21
type_of_loss =      1
!---------------------------Kernel / -------------------------!
write_kernel_matrix = .false. 
kernel_type =      4
kernel_rff =      1
np_kernel_ref =    200
np_kernel_full =    800
np_omega =      4
kernel_dump =      1
length_kse =  0.100000000000000D+01
sigma_kse =  0.100000000000000D+01
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
selection_type =     3
seed =  14792435
weighted = .false. 
weighted_auto = .false. 
weighted_3ch = .false. 
fix_no_of_elements =     1
iread_energy =     2
chemical_elements = " Fe"
chemical_elements_invisible = ""
weight_per_element = "1.d0"
weight_per_element_3ch = "1.d0"
ref_energy_per_element = "0.d0"
n_pca =     3
classes_for_mcd = " 01 02"
sign_stress =  0.100000000000000D+01
sign_stress_big_box =  0.100000000000000D+01
zbl_potential = .false. 
r1_zbl =  0.100000000000000D+01
r2_zbl =  0.220000000000000D+01
force_variance = .false. 
energy_variance = .false. 
stress_variance = .false. 
!------------------------Descriptors------------------------!
descriptor_type =      9
write_desc = .true. 
lmask = .false. 
mask_file = .false. 
write_desc_dump = .false. 
read_desc_dump = .false. 
desc_forces = .false. 
val_desc_max =  0.100000000000000D+01
r_cut =  0.470000000000000D+01
r_cut_in = -0.120000000000000D+01
r_cut_width_in =  0.400000000000000D+00
r_cut_width =  0.500000000000000D+00
type_fcut =      2
!---G2----
n_g2_eta =      3
n_g2_rs =      1
eta_max_g2 =  0.200000000000000D+03
eta_min_g2 =  0.200000000000000D+00
rs_max_g2 =  0.200000000000000D+03
rs_min_g2 =  0.200000000000000D+00
!---G3----
n_g3_eta =      3
n_g3_zeta =      2
n_g3_lambda =      2
!---BehPar----
strict_behler = .false. 
!---AFS----
afs_type =      1
n_rbf =      4
n_rbf_afs =      4
n_rbf_so3 =      4
n_cheb =      5
!---MTP----
mtp_poly_min =      2
mtp_poly_max =      5
!---SO3/BS03----
l_max =      4
radial_pow_so4 =      1
!---SO4/BS04----
j_max =    4.0000000000
inv_r0_input =  0.993633802276324D+00
lbso4_diag = .false. 
!---SOAP----
alpha_soap =    2.0000000000
n_soap =      2
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
l_body_order(3) = .true. 
l_body_order(4) = .true. 
l_body_order(5) = .false. 
bond_dist_transform =      3
bond_beta =  0.200000000000000D+01
bond_dist_ann =  0.100000000000000D+01
!----MiLaDy----
rmat_dim =     20
img_num_ch =      1
power_line =  0.200000000000000D+01
power_coeff_renorm =  0.333000000000000D+00
img_weighted = .false. 
ws1 = "1.d0"
ws2 = "1.d0"
ws3 = "1.d0"
ws4 = "1.d0"
!----mag----
magnetic_sld_j_dim =     10
magnetic_sld_s2_dim =     10
magnetic_sld_s4_dim =     10
!---FOURIER NBODY----
dim_fourier_nbody = " 0 0 2000 0 0"
length_fourier_nbody = " 0.7d0 0.7d0 0.7d0 0.7d0 0.7d0"
sigma_fourier_nbody = " 1.d0 1.d0 1.d0 1.d0 1.d0"
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
optimize_weights_Le = .false. 
optimize_ga_population =     40
class_no_optimize_weights = ""
max_iter_optimize_weights =     40
factor_energy_error =  0.100000000000000D+01
factor_force_error =  0.100000000000000D+01
factor_stress_error =  0.100000000000000D+01
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
n_kcross =         0
dim_data =       100
nd_fingerprint =         3
max_data =      0
n_kcross =      0
& end 
