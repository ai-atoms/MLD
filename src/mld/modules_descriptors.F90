!
!module descriptors_local
!  !use ondm_gen_com_m, only : imm
!  !use ml_in_ndm_module, only: imm_neigh
!
!  double precision, dimension(:, :), allocatable     :: local_g2
!  double precision, dimension(:, :, :, :), allocatable     :: local_g2_deriv
!
!  double precision, dimension(:, :), allocatable     :: local_g3
!  double precision, dimension(:, :, :, :), allocatable     :: local_g3_deriv
!
!  double precision, dimension(:, :), allocatable     :: local_pow_so3
!  double precision, dimension(:, :, :, :), allocatable     :: local_pow_so3_deriv
!
!  double precision, dimension(:, :), allocatable     :: local_bispectrum_so3
!  double precision, dimension(:, :, :, :), allocatable     :: local_bispectrum_so3_deriv
!
!  double precision, dimension(:, :), allocatable     :: local_pow_so4
!  double precision, dimension(:, :, :, :), allocatable     :: local_pow_so4_deriv
!
!  double complex, dimension(:, :), allocatable :: local_bispectrum_so4
!  double complex, dimension(:, :, :, :), allocatable :: local_bispectrum_so4_deriv
!
!  double precision, dimension(:, :), allocatable     :: local_soap
!  double precision, dimension(:, :, :, :), allocatable     :: local_soap_deriv
!
!  double precision, dimension(:, :), allocatable     :: local_mtp
!  double precision, dimension(:, :, :, :), allocatable     :: local_mtp_deriv
!
!  double precision, dimension(:, :), allocatable     :: local_dmilady
!  double precision, dimension(:, :, :, :), allocatable     :: local_dmilady_deriv
!
!  double precision, dimension(:, :), allocatable     :: local_afs
!  double precision, dimension(:, :, :, :), allocatable     :: local_afs_deriv
!
!  !imm size
!  integer, dimension(:), allocatable     :: d_n_neigh, l_d_n_neigh
!  !imm, imm_neigh size
!  integer, dimension(:, :), allocatable  :: l_d_kind_neigh, d_kind_neigh
!
!end module descriptors_local
