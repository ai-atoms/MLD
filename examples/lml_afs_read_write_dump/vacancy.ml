&input_ml
debug=.false.

!#####################################################################
!                         DEFINE THE ML TASK
!#####################################################################

ml_type=0            ! -1 - only compute descriptors without any fit
                     !  0 - functions SNAP etc,
                     !  1 - KRR
                     !  2 - GAP


mld_order=1         ! 1 - linear , 2 quadratic ; 11 - n-linear

! -- DIFFERENT NUMERICAL WAYS TO PERFORM THE FIT ---

mld_fit_type=3      ! 0  - home made BH ;
                     ! 10 - home made with regularization ;
                     ! 1  - lapack QR ;
                     ! 2  - lapack + constraints QR based snap-class_constraint should be written  ;
                     ! 3  - lapack full ortho decomposition with rank estimation. Norm minimization.
                     ! 4  - lapack full SVD  with rank estimation. Norm minimization.

!--- WRITING AND READIND DESCRIPTORS ---

write_desc=.false.
desc_forces=.true.

!--Dump options--
write_desc_dump = .false.
read_desc_dump = .true.

!--- HOW TO SELECT THE FILES FROM DB ---

selection_type = 1   ! 1 - selects first N_train elements of the database;
                     ! 2 - selects last  N_train elements of the database;
                     ! 3 - selects randomly N_train elements of the database.

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

!--- AFS ---
afs_type=1          ! set 1 for standard AFS; set 2 for modified AFS with tensorial product in the radial channels.
n_rbf=12            ! The number of radial channels
n_cheb=8            ! The number of Chebyshev polynomials


! For the settings of other descriptors, see the manual



&end
