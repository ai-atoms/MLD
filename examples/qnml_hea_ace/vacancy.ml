&input_ml

!--------------------------------------------------------------------\
!                           ML Mode                                  |
!--------------------------------------------------------------------/
ml_type=0

!---------------ML model--------------------------!
mld_order=2             ! 1 - linear , 2 quadratic , 11 - n-linear, 3 PolyC, 7 - kernel, 11 - nlinear
mld_type_quadratic=1    ! 1 default +LML precondition, 2 - full quadratic, 3 - bi-linear, 4 - zx , 5 - ZX

desc_forces=.true.

!--------------------------------------------------------------------\
!                           ML Train                                 |
!--------------------------------------------------------------------/
mld_fit_type=4

!--- for the systems with several elements ------
weighted=.true. ! if true : weighted descriptors 
weighted_3ch=.false.
fix_no_of_elements=4
chemical_elements=" Ta Ti V  W"
weight_per_element="0.8 0.9 1.0 1.1"  ! for numerical stability, set the values close to 1.0

!--> cut-off of the descriptor ...
r_cut=4.7d0
r_cut_width=0.5d0 
r_cut_in=1.2d0 
r_cut_width_in=0.4d0 
type_fcut=3

!---> Type of  descriptors 
descriptor_type=300 ! 300 - ACE/kACE

ace_numax=1
ace_gencg=2  ! 1 - is DRAFT redundant version ; 2  - is SVD Dusson-Ortner version
ace_chem = 1 ! 0 - incomplet ; 1 - standard ; 2 - TS
ace_radial_chem = 1 ! 1 - Ralf (standard ACE) ; 3 - HSVD (kACE) ; 5 - HSVD with random projection

l_ace_order(1)=.true.
l_ace_order(2)=.false.
l_ace_order(3)=.false.
l_ace_order(4)=.false.
l_ace_order(5)=.false.
l_ace_order(6)=.false.
ace_nmax_list="2 1 1 1 1 1"
ace_lmax_list="0 2 2 2 1 1"
ace_lambda_list="3.0 3.0 3.0 3.0 3.0 3.0"
ace_radial_poly = 2 ! 1 - powPftouny ; 2 - expPaftouny ; 3 - simpBessel 

&end
