&input_ml

debug=.false.
debug_time=.false.
!--------------------------------------------------------------------\
!                           ML Mode                                  |
!--------------------------------------------------------------------/


ml_type=0            ! -1 - just compute descriptors ...
                     !  0 - functions SNAP etc,
                     !  1 - Kernel RR
                     !  2 - GAP, not yet implemented
                     !  3 - n - linear
                     ! -2 - analyse the data or dump kernel.
!-------------------------------------------------:


!--------------------------------------------------------------------\
!                           ML Model                                 |
!--------------------------------------------------------------------/

!---------------ML model--------------------------!
snap_order=1             ! 1 - linear , 2 quadratic , 11 - n-linear, 3 PolyC, 7 - kernel, 11 - nlinear
snap_type_quadratic=1    ! 1 default +LML precondition, 2 - full quadratic, 3 - bi-linear
!this is for n-linear
order_nlinear = 2 ! is the order of n-linear form. Active if snap_order=11. For the moment is implemented only the case 2.
!this is for PolyC
polyc_n_poly = 2    ! max Poly: 1(lieanr) to 3
polyc_n_hermite = 4 ! Max Hermite degree: 1(identity) to 4

!--------------------------------------------------------------------\
!                      DB visualization                              |
!--------------------------------------------------------------------/
n_pca=3
classes_for_mcd="07 08 09 13 11"

!--------------------------------------------------------------------\
!                           Database parameters                      |
!--------------------------------------------------------------------/
!build_subdata=.false.
selection_type=1 ! 1 - selects first    "ns" elements of the database
                 ! 2 - selects last     "ns" elements of the database
                 ! 3 - selects randomly "ns" subsets of "kelem" elements of the database
seed=37          ! seed for random generator
db_file="db_model.in"
db_path="DB/"

!some compatibilities for stress betweene databases ....
!sign_stress=0.3333333d0
! W GABOR is -
!sign_stress=-1.d0
! JUlien is 1
sign_stress=1.d0
sign_stress_big_box=1.d0


!--------------------------------------------------------------------\
!                      KERNEL                                        |
!--------------------------------------------------------------------/
!write_kernel_matrix=.true.
power_mcd = 0.10d0
np_kernel_ref = 1000
np_kernel_full = 2000
kernel_dump=3 ! 1 - dump by classes_for_kernel.
              ! 2 - dump by MCD distances. !NOT YET
              ! 3 - dump_by MAHALANOBIS
! 1 Gaussian ; 4 Polynomial ; MCD like
kernel_type=6   ! 1 - squared exponential OK
                ! 2 - orstein uhlenbeck
                ! 3 - matern class
                ! 4 - polynomial     OK
                ! 5 - MCD            OK
                ! 6 - MAHA           OK
sigma_kernel=0.0
length_kernel=0.005
length_kernel=0.707106
kernel_power=0.5

!If you want to read the desgn matrix.
write_design_matrix = .false.
!-------------------------------------------------:



!--------------------------------------------------------------------\
!                           ML Train                                 |
!--------------------------------------------------------------------/

! if only training:
train_only=.false.

!--------------train options----------------------!
!snap_fit_type
!0  - home made BH ;
!10 - home made with regularization ;
!1 - lapack QR with norm minimization ;
!2 - lapack + constraints QR based snap-class_constraint should be written  ;
!3 - lapack full ortho decomposition with rank estimation. Norm minimization.
!4 - lapack full SVD  with rank estimation. Norm minimization.
snap_fit_type=1
svd_rcond=1.d-14
snap_regularization_type=0  ! 1 - grid regularization ; 0 - none
snap_class_constraints="06"


!--------------regularization---------------------!
!regularization ....
!with j=3.5 was 1.d-4 (default value)
! with 1.d-8 and SVD good results for QNML
lambda_krr=1.d-10
!lambda_krr=0.3353550443E+00
min_lambda_krr= 1.d-12
max_lambda_krr= 1.d+01
n_values_lambda_krr=20

!weights optimization using genetic algorithms ...
optimize_ga_population=80
optimize_weights=.false.
class_no_optimize_weights="02"
max_iter_optimize_weights=20
optimize_weights_L1=.false.
optimize_weights_L2=.false.
optimize_weights_Le=.true.
factor_energy_error=1.d0
factor_force_error=1.d0
factor_stress_error=1.d0
!--------------------------------------------------!





!--------------------------------------------------------------------\
!                        Descriptors parameters                      |
!--------------------------------------------------------------------/

weighted=.true. ! if true : weighted descriptors and number of type of descriptors n2_typ,n3_typ = 1
weighted_auto=.true.
weighted_3ch=.true.
fix_no_of_elements=4
chemical_elements=" Ta Ti V W  "
!chemical_elements="  "
weight_per_element="0.1E+01  0.1E+01  0.1E+01  0.1E+01 "
weight_per_element_3ch="0.1E+01  0.1E+01  0.1E+01  0.1E+01 "
!weight_per_element="0.1249500920E+01  0.6426548602E+00  0.6629758219E+00  0.1259446403E+01 "
!weight_per_element_3ch="0.1049500920E+01  0.2426548602E+00  0.1629758219E+00  0.2759446403E+01 "
!weight_per_element=""

!--> store by writting descriptors ....
write_desc=.false.
!w write_desc_dump=.false.
!w read_desc_dump=.false.
!--> nodescriptors for forces ... unusefull is no forces / stress are fitted.
desc_forces=.true.

!--> cut-off of the descriptor ...
r_cut=5.0d0



!val_desc_max=1.d0

!---> The descriptord type
descriptor_type=4   ! 1 - g2
                    ! 2 - g3
                    ! 3 - Behler
                    ! 4 - Angular fourier series (afs)
                    ! 5 - soap (debug)
                    ! 6 - power spectrum SO3
                    ! 7 - bispectrum SO3
                    ! 8 - power spectrum SO4
                    ! 9 - bispectrum SO4
                    ! 14 - G2  + AFS
                    ! 18 - G2  + pSO4
                    ! 19 - G2  + bispectrum SO4
                    ! 30 - Magnetic_SLD energy only
                    ! 34 - Magnetic_SLD_AFS  energy only
                    ! 100 - MTP
                    ! 200 - body PiP




!---> desc_milady TM
rmat_dim=50


!--->  set-up MTP
mtp_poly_min=2
mtp_poly_max=7

!--->  set-up G2
n_g2_eta=15          !g2
n_g2_rs=1            !g2
eta_max_g2=2.0

!--->  set-up G3
n_g3_eta=7          ! g3
n_g3_zeta=3         ! g3
n_g3_lambda=2       ! g3

strict_behler=.true.

body_D_max(2)=14
body_D_max(3)=12
body_D_max(4)=4
body_D_max(5)=2
l_body_order(2)=.true.
l_body_order(3)=.true.
l_body_order(4)=.true.
l_body_order(5)=.false.
! 1 - dist, 2 - exp 3 - inverse
bond_dist_transform=3
bond_beta=2.0
bond_dist_ann=1.0

!---> AFS
afs_type=1
! 1 is the standard from PRB 2003 . With two there is tensorial product in the radial channels.
n_rbf=8            ! is used by afs, pow_so3, bispectrum_so3
n_cheb=8           ! afs

!---> BODY




!---> set up SO4 and bi-SO4
j_max=4.0           ! pow_so4, bispectrum_so4
lbso4_diag=.false.  !lbso4_diag=.true.  ... only diagonal (Gabor like)
                    !lbso4_diag=.false. ... only diagonal (non diagonal -SNAP like)


!---> SOAP only ...
!soap=.false.
alpha_soap=3.0d0        ! soap
n_soap=7                ! soap radial ...
lsoap_diag=.false.      !
lsoap_norm=.true.       !
lsoap_lnorm=.true.      !
r_cut_width_soap=0.5d0  !
!soap_fcut=.true.       ! soap


! l_max is used by pow so3 and bi so3 and SOAP
!bispectrum_so3: doesn't forget that bso3 is defined by l_max and n_rbf.
l_max=4
                        ! lbso3_diag=.false. ... only diagonal (non diagonal -SNAP like)



!---------------------------------------------------------------------!
!              .............. LBFGS setting .............             !
!---------------------------------------------------------------------!
lbfgs_xtol=1.d-20
lbfgs_eps=1.d-2
lbfgs_m_hess=100
lbfgs_max_steps=10000
lbfgs_print(1)=100
lbfgs_print(2)=0
lbfgs_gtol=9.d-01


length_kse=1.d0
sigma_kse=100.1d0

iread_ml=0
isave_ml=0

nd_fingerprint=3

!----------------------------------------------------------------------------------

&end
