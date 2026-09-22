&input_ml
debug=.false.


ml_type=0            ! -1 - only compute descriptors ...
                     ! 0 - functions SNAP etc,
                     ! 1 KRR
                     ! 2 GAP


mld_order=2             ! 1 - linear , 2 quadratic ; 11 - n-linear
mld_type_quadratic=1     ! 1 default +LML precondition, 2 - full quadratic, 3 - bi-linear


desc_forces=.true.

!sign_stress=-1.d0
!sign_stress_big_box=1.d0


!--------------train options----------------------!


mld_fit_type=4  ! 0  - home made BH ;
                 ! 10 - home made with regularization ;
                 ! 1  - lapack QR ;
                 ! 2  - lapack + constraints QR based snap-class_constraint should be written  ;
                 ! 3  - lapack full ortho decomposition with rank estimation. Norm minimization.
                 ! 4  - lapack full SVD  with rank estimation. Norm minimization.





!--- for the systems with several elements ------

weighted=.true. ! if true : weighted descriptors 
weighted_3ch=.false.
fix_no_of_elements=4
chemical_elements=" Ta Ti V  W"
weight_per_element="0.8 0.9 1.0 1.1"  ! for numerical stability, set the values close to 1.0


!---------------------------------------------------------------------
! Descriptor parameters
!----------------------------------------------------------------------



!--> cut-off of the descriptor ...
r_cut=4.7d0


!---> Type of  descriptors 
descriptor_type=9   ! 1 - g2
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



!---> set-up SO4 and bi-SO4
j_max=2.5           ! j=3.5 is a reasonable choice; 


&end
