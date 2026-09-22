&input_ml
debug=.false.

!#####################################################################
!                         DEFINE THE ML TASK
!#####################################################################

ml_type=1            ! -1 - only compute descriptors without any fit
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

!---> BSO4 
j_max = 4

!--------------------------------------------------------------------\
!                      KERNEL                                        |
!--------------------------------------------------------------------/
! ----- for sparsification if needed
write_kernel_matrix=.false.
power_mcd = 0.05d0
kernel_dump=3 ! 1 - dump by classes_for_kernel.
              ! 2 - dump by MCD distances. !NOT YET
              ! 3 - dump_by MAHALANOBIS
! ----- kernel details
np_kernel_ref = 2000
np_kernel_full = 2000

! 1 Gaussian ; 4 Polynomial ; MCD like
kernel_type=7   ! 1 - squared exponential OK se
                ! 2 - orstein uhlenbeck
                ! 3 - matern class
                ! 4 - polynomial     OK  po
                ! 5 - MCD            not OK
                ! 6 - MAHA           OK
                ! 7 - random kernel  OK  random
sigma_kernel=0.02
!length_kernel=0.005
!length_kernel=0.707106


! For the settings of other descriptors, see the manual



&end
