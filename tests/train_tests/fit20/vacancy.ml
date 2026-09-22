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
descriptor_type=202

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
