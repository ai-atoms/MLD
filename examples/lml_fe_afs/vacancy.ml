&input_ml
debug=.false.

!#####################################################################
!                         DEFINE THE ML TASK
!#####################################################################

ml_type=0            ! -1 - only compute descriptors without any fit
                     !  0 - functions SNAP etc,
                     !  1 - KRR
                     !  2 - GAP


mld_order=2                 ! 1 - linear , 2 quadratic ; 11 - n-linear


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

!---> AFS
afs_type=1         ! set 1 for standard AFS; set 2 for modified AFS with tensorial product in the radial channels.
n_rbf=4            ! The number of radial channels
n_cheb=4           ! The number of Chebyshev polynomials


! For the settings of other descriptors, see the manual

&end
