&input_ml
debug=.false.
scalapack_driver=.true.
nbr_predefined=64
nbc_predefined=64
!--------------------------------------------------------------------\
!                           ML Model                                  |
!--------------------------------------------------------------------/
ml_type=1            !set 1 to perform the fit using kernel
snap_order=1         !set 1 for linear regression
snap_fit_type=1      !scalapack full QR  without rank estimation
lambda_krr=1.d-8     ! regularization (default)
desc_forces=.true.
!--------------------------------------------------------------------\
!                          Your system                                |
!--------------------------------------------------------------------/

fix_no_of_elements=1
chemical_elements=" W "
weight_per_element="1.d0"

!--------------------------------------------------------------------\
!                        Descriptors parameters                      |
!--------------------------------------------------------------------/

!--> cut-off of the descriptor ...
r_cut=5.0d0
descriptor_type=9   ! 9 - bispectrum SO4
j_max=4.0           ! bispectrum_so4
lbso4_diag=.false.  !lbso4_diag=.true.  ... only diagonal (Gabor like)
                    !lbso4_diag=.false. ... only diagonal (non diagonal -SNAP like)



!--------------------------------------------------------------------\
!                      DB kernel setting                             |
!--------------------------------------------------------------------/
write_kernel_matrix=.false.          ! write or no the selected kernel. 
n_pca=3                             ! used for vizualization (default) 
classes_for_mcd="07 08 13"          ! classes used for the MCD reference
power_mcd = 0.05d0                  ! exponenet used for the grid of selction of kernel
np_kernel_ref = 1000                ! number of proposed points in the MCD class 
np_kernel_full = 1000               ! number of points outside the MCD class
kernel_dump=3                       ! set 3 - dump_by MCD/MAHALANOBIS


!--------------------------------------------------------------------\
!                      kernel hyper-param                            |
!--------------------------------------------------------------------/
kernel_type=4      ! 4 - polynomial     OK  po
sigma_kernel=0.0   ! set 0 for polynomial in milady
length_kernel=0.05 ! no known value test many.
kernel_power=4     ! 4 - good compromise. order of polynomial kernel 1 linear, 2 quadratic (default) 

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


&end
