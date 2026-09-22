&input_ml
debug=.false.


ml_type=0            ! -1 - just compute descriptors ...
                     ! 0 - functions SNAP etc,
                     ! 1 KRR
                     ! 2 GAP, not yet implemented




!snap_fit_type
!0  - home made BH ;
!10 - home made with regularization ;
!1 - lapack QR ;
!2 - lapack + constraints QR based snap-class_constraint should be written  ;
!3 - lapack full ortho decomposition with rank estimation. Norm minimization.
!4 - lapack full SVD  with rank estimation. Norm minimization.
snap_fit_type=4
svd_rcond=-1
snap_regularization_type=0   ! 1 - grid regularization ; 0 - none
snap_class_constraints="06"


lambda_krr=1.d-8





weighted=.false. ! if true : weighted descriptors and number of type of descriptors n2_typ,n3_typ = 1
fix_no_of_elements=1
chemical_elements=" W "



!--------------------------------------------------!

!--> cut-off of the descriptor ...
r_cut=5.00d0
r_cut_width=0.5d0 
r_cut_in=1.2d0 
r_cut_width_in=0.4d0 
type_fcut=3 
!val_desc_max=1.d0

! ZBL ---------------------
zbl_potential = .false.
r1_zbl=1.0d0
r2_zbl=2.2d0
activate_k2b = .false.
sigma_2b=0.3d0
delta_2b=1.d0
np_radial_2b=30



!---------------------------------------------------------------------
! Descriptor parameters
!----------------------------------------------------------------------

!--> store by ritting descriptors .... 
write_desc=.false.
!--> nodescriptors for forces ... unusefull is no forces / stress are fitted. 
desc_forces=.true. 

!--> cut-off of the descriptor ...
r_cut=5.0d0


!---> The descriptord type 
descriptor_type=300  ! ACE 
ace_numax = 6
ace_chem = 1 ! 0 - incomplet ; 1 - standard 
ace_radial_chem = 3 ! 1 - Ralf ; 2 - BLOCK HSVD  ; 3 - HSVD ; 4 - CHEM MAPS 
ace_gencg= 1  ! 1 - is DRAFT redundant version ; 2  - is SVD Dusson-Ortner version
l_ace_order(1)=.true.
l_ace_order(2)=.true.
l_ace_order(3)=.true.
l_ace_order(4)=.true.
l_ace_order(5)=.true.
l_ace_order(6)=.false.
!REF1 is for ace_numax=3
ace_kmax_list="8  8 8 8 5  1"
ace_nmax_list="14 4 3 2 1 1 "
ace_lmax_list="0  3 2 2 1 1"
ace_lambda_list="3.0 3.0 3.0 3.0 3.0 3.0"
ace_radial_poly = 2 ! 1 - powPftouny ; 2 - expPaftouny ; 3 - simpBessel 
l_ace_set_rcut = .false. 
ace_rcut_in_list = " 1.2d0 1.2d0 1.2d0 1.2d0 1.2d0 1.2d0"  
ace_rcut_out_list = " 5.d0 5.d0 5.d0 5.d0 5.d0 5.d0"  
ace_rcut_width_in_list = " 0.4d0 0.4d0 0.4d0 0.4d0 0.4d0 0.4d0"
ace_rcut_width_out_list = " 0.5d0 0.5d0 0.5d0 0.5d0 0.5d0 0.5d0"

zetaace_order = 1
! here is the body nu=1... 
dim_delta_zetaace(1) = " 1.d0 1.d0 1.d0 1.d0 1.d0 1.d0 "  
! nu =2 
dim_delta_zetaace(2) = " 1.d0 1.0d0 1.0d0 1.0d0 1.d0 1.d0 "  
! nu = 3 
dim_delta_zetaace(3) = " 1.d0 1.0d 1.0d0 1.0d0 1.d0 1.d0 "  
! nu = 4 
dim_delta_zetaace(4) = " 1.d0 1.0d 1.0d0 1.0d0 1.d0 1.d0 "  


body_D_max(2)=20
body_D_max(3)=9
body_D_max(4)=7
body_D_max(5)=4

l_body_order(1)=.false.
l_body_order(2)=.true.
l_body_order(3)=.true.
l_body_order(4)=.false.
l_body_order(5)=.false.
! 1 - dist, 2 - exp 3 - inverse
bond_dist_transform=3
bond_beta=2.0
bond_dist_ann=2.0

!find_best_length = .true.
r_cut_ft2b = 5.3d0
r_cut_width_ft2b = 1.5d0
r_cut_ft3b = 5.3d0
r_cut_width_ft3b = 1.5d0
r_cut_ft4b = 4.6d0
r_cut_width_ft4b = 1.3d0
r_cut_ft5b = 3.3d0
r_cut_width_ft5b = 0.9d0


!random
sigma_kernel = 1.0d0 
length_kernel=0.5d0
krff_type=1  ! 1 Gaussian;  2 Cauchy;  3 Laplace 
dim_fourier_nbody="0 80 1000 2000 4  "
!kernel_power=0.02




! build subdatabase
!----------------------------------------------------------------------------------
selection_type=1 ! 1 - selects first    "ns" elements of the database
                 ! 2 - selects last     "ns" elements of the database
                 ! 3 - selects randomly "ns" subsets of "kelem" elements of the database
seed=37          ! seed for random generator
db_file="db_model.in"
db_path="DB/"


!----------------------------------------------------------------------------------

&end
