&input_ml
debug=.false.

!#####################################################################
!                         DEFINE THE ML TASK
!#####################################################################

ml_type=0            ! -1 - only compute descriptors without any fit
                     !  0 - functions SNAP etc,
                     !  1 - KRR
                     !  2 - GAP


mld_order=1          ! 1 - linear , 2 quadratic ; 11 - n-linear


write_desc=.false.
desc_forces=.true.



!#####################################################################
!                    DEFINE YOUR ATOMIC SYSTEMS
!#####################################################################

weighted=.false.
chemical_elements=" Fe "

!#####################################################################
!                         DEFINE YOUR DESCRIPTORS
!#####################################################################

r_cut=4.7d0

descriptor_type=9   ! 9 - bispectrum SO4


! BSO4
j_max=4             ! 



!#####################################################################
!                         DEFINE YOUR K2B
!#####################################################################
activate_k2b = .true.
sigma_2b = 0.3d0          ! controls the smoothness of the interaction
delta_2b = 1.d0           ! controls the magnitude (in eV) of a typical 2-body interaction
np_radial_2b = 40         ! number of radial points in the 2-body grid
r_cut_in = 1.5d0          !  where the many-body interaction is shut dow (in A)n, and a smooth link with ZBL is prepared. Below this distance, there is no fitting.
r_cut_width_in = 0.3d0    ! controls the smoothness of the inner cutoff
type_fcut = 3             ! the type of fcut function (there are several options).  Iforgot if the best is 2 or 3 :) 


!#####################################################################
!                         ACTIVATE YOUR ZBL
!#####################################################################

zbl_potential = .true.    ! activate or not
r1_zbl = 1.0d0             ! for r < r1_zbl, ONLY ZBL is used in the radial part of your ML
r2_zbl = 2.2d0             ! internal parameter that controls the smoothness of ZBL



&end
