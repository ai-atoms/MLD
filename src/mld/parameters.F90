! HND XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
! HND X
! HND X   MiLaDy (Machine Learning Dynamics) 
! HND X   
! HND X
! HND X   Milady was written and designed by Mihai-Cosmin Marinica, Alexandra M. Goryaeva 
! HND X   Copyright 2015-2025.
! HND X
! HND X   Portions of MiLady were written by Wesley UnnToc, Clovis Lapointe, Anruo Zhong, 
! HND X   Jacopo Baima, Anida Khizar, Christian van Wambeke 
! HND X
! HND X   MiLaDy is published and distributed under the
! HND X      Academic Software License v1.0 (ASL)
! HND X
! HND X   MiLaDy is distributed at https://github.com/ai-atoms/milady in the hope that it 
! HND X   will be useful for non-commercial academic research, but WITHOUT ANY WARRANTY;  
! HND X   without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
! HND X   PURPOSE. See the ASL for more details.
! HND X
! HND X   You should have received a copy of the ASL along with this program
! HND X   (e.g. in a LICENSE.md file); if not, you can write to the original licensors,
! HND X   Mihai-Cosmin Marinica and Alexandra M. Goryaeva. The ASL is also published at
! HND X   http://github.com/ai-atoms/milady/ASL
! HND X
! HND X   When using this software, please cite the following reference:
! HND X
! HND X   A.M. Goryaeva et al. Comput. Mater. Sci. 166, 200-209 (2019)
! HND X   A.M. Goryaeva et al. Phys. Rev. Mater. 5: 103803 (2021)
! HND X
! HND XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

#include "../MLD_MACROS.INC"

subroutine train_snap_get_parameters_linear_scalapack()

  use ml_in_ndm_module, only: debug, rangml, &
                              mld_fit_type, &
                              fit_home_hb, fit_lapack_qr,  &
                              fit_lapack_ortho, fit_lapack_svd, fit_als, fit_online_svd, svd_rcond, &
                              descriptor_type, descriptor_body, descriptor_ace, &
                              lambda_krr
  use module_input_als_fit, only: als_nsteps, als_tol, als_ridge_k, als_rho, &
                              als_alpha_method, als_precond_type, als_precond_flat, &
                              als_nnls_alpha, als_nnls_mode, &
                              als_nu_max, als_block_partition, als_alpha, &
                              set_als_block_partition

  use module_fit_ScaMatrix, only: prepare_sca_phi, scalapack_lsystem_by_homeLU, & 
                                  scalapack_lsystem_by_QR, & 
                                  scalapack_square_lsystem_by_cholesky, &
                                  scalapack_lsystem_by_SVD
  use module_para_als, only: para_als_ridge_solve

  use module_ml_scalapack, only: context, nbr_Amat, dimr_sca_Amat, &
                                 dimr_sca_phi, dimc_sca_phi, sca_phi, desc_sca_phi, nbr_phi, &
                                 nbc_phi, sca_ymat_qr_svd, desc_sca_ymat_qr_svd, dimr_sca_ymat_qr_svd, &
                                 dimc_sca_ymat_qr_svd, &
                                 sca_w_params, desc_sca_w_params, dimr_sca_w_params, dimc_sca_w_params, &
                                 l_dimc_sca_w_params, l_dimr_sca_w_params, nbr_w_params, nbc_w_params,  &
                                sca_ymat, desc_sca_ymat, dimr_sca_ymat,& 
                                myrow, mycol, nprow, npcol, iam, &
                                sca_Amat, desc_sca_Amat, dimc_sca_Amat, nbc_Amat, &
                                l_dimr_sca_Amat, l_dimc_sca_Amat
  use module_scalapack_tools, only:   allocation_scalapack_matrix

  use snap, only: w_params, ymat, weights_snap
  use mld_logger

  implicit none
  integer :: rank_sca
  integer :: ii
  _NAMECURRENT_("train_snap_get_parameters_linear_scalapack")



  _MLD_BEGIN_

  ! ToScaRemove
  if (debug) then 
  if (rangml==0) write(6,*) 'Before prepare_sca_phi'
  end if 
  ! end ToScaRemove
  ! sca_w_params sca version for w_params
  !         |      ...   |
  !         |      ...   |
  !       = |  f(D) x 1  |  f(D) is line numbers of Amat
  !         |      ...   |
  !         |      ...   |
  nbr_w_params = nbr_Amat 
  nbc_w_params = 1
  dimr_sca_w_params = dimr_sca_Amat
  dimc_sca_w_params = 1
  call allocation_scalapack_matrix (" sca_w_params ", sca_w_params, desc_sca_w_params, &
                                  dimr_sca_w_params, dimc_sca_w_params, &
                                  l_dimr_sca_w_params, l_dimc_sca_w_params, nbr_w_params, nbc_w_params, & 
                                  myrow, mycol, nprow, npcol, context, iam  )

  ! set sca_ymat from ymat ... 
  do ii = 1, dimr_sca_ymat
    call pdelset(sca_ymat, ii, 1, desc_sca_ymat, ymat(ii,1)) 
  end do 

  call prepare_sca_phi
  call blacs_barrier(context, 'A')

  ! ToScaRemove
  if (debug) then 
     if (rangml==0) write(6,*) 'Out from  prepare_sca_phi'
  end if    
  ! end ToScaRemove



  select case(mld_fit_type)

    !this is equivalent to HOME - only square matrix - Working tested. 
    case(fit_home_hb)

     call scalapack_lsystem_by_homeLU(dimr_sca_phi, dimc_sca_phi, sca_phi, desc_sca_phi, &
                                             nbr_phi, nbc_phi, &
                                             dimr_sca_ymat_qr_svd, dimc_sca_ymat_qr_svd, &
                                             sca_ymat_qr_svd, desc_sca_ymat_qr_svd, w_params)

    !this is QR -  general MxN matrix but full rank - Working tested.
    case (fit_lapack_qr)
      call scalapack_lsystem_by_QR(dimr_sca_phi, dimc_sca_phi, sca_phi, desc_sca_phi, & 
                                  nbr_phi, nbc_phi, &
                                  dimr_sca_ymat_qr_svd, dimc_sca_ymat_qr_svd, &
                                  sca_ymat_qr_svd, desc_sca_ymat_qr_svd, w_params)

    !this is Cholesky - only square matrix positive definite. 
    case (fit_lapack_ortho)
      call scalapack_square_lsystem_by_cholesky (dimr_sca_phi, dimc_sca_phi, sca_phi, desc_sca_phi, & 
                                  nbr_phi, nbc_phi, &
                                 dimr_sca_ymat_qr_svd, dimc_sca_ymat_qr_svd, &
                                 sca_ymat_qr_svd, desc_sca_ymat_qr_svd, w_params)
  
    !this is home made SVD for LSS problems, general matrix
    case(fit_lapack_svd)
      call scalapack_lsystem_by_SVD (dimr_sca_phi, dimc_sca_phi, sca_phi, desc_sca_phi, & 
                           nbr_phi, nbc_phi, &
                           dimr_sca_ymat_qr_svd, dimc_sca_ymat_qr_svd, &
                           sca_ymat_qr_svd, desc_sca_ymat_qr_svd, w_params, &
                           svd_rcond, rank_sca, context) 

    !mld_fit_type=5 ALS-Ridge block preconditioning (ScaLAPACK path)
    case(fit_als)
      ! ALS-Ridge: works directly on the distributed design matrix sca_Amat
      ! bypasses the phi = A W A^T normal equation construction.
      ! Apply sqrt(weights_snap) to sca_Amat columns and sca_ymat rows in-place
      ! so the ALS solver minimises sum_j w_j ||A_j w - y_j||^2.
      ! sca_Amat is D x M distributed by columns (block nbc_Amat, mycol).
      ! sca_ymat is M x 1 distributed by rows  (block nbc_Amat, mycol) - same j-index.
      block
        integer :: ilocal_w, jlocal_w, jglobal_w
        integer :: INDXL2G
        real(kind=kind(1.d0)) :: sw, val_w
        do jlocal_w = 1, l_dimc_sca_Amat
          jglobal_w = INDXL2G(jlocal_w, nbc_Amat, mycol, 0, npcol)
          sw = dsqrt(weights_snap(jglobal_w))
          do ilocal_w = 1, l_dimr_sca_Amat
            sca_Amat(ilocal_w, jlocal_w) = sca_Amat(ilocal_w, jlocal_w) * sw
          end do
          call pdelget('A', ' ', val_w, sca_ymat, jglobal_w, 1, desc_sca_ymat)
          call pdelset(sca_ymat, jglobal_w, 1, desc_sca_ymat, val_w * sw)
        end do
      end block
      ! Resolve the block partition from descriptor structure (or keep user-supplied)
      call resolve_als_partition_from_descriptors(dimr_sca_Amat)
      block
        real(kind=kind(1.d0)), allocatable :: alpha_als(:)
        integer :: als_info, inu
        if (myrow == 0 .and. mycol == 0) then
          do inu = 1, als_nu_max
            write(6, '("ALS-Ridge ScaLAPACK: block ", i4, " rows [", i8, ",", i8, "]")') &
              inu, als_block_partition(2*inu-1), als_block_partition(2*inu)
          end do
        end if
        allocate(alpha_als(als_nu_max))
        call para_als_ridge_solve(sca_Amat, desc_sca_Amat, &
                                  dimr_sca_Amat, dimc_sca_Amat, &
                                  nbr_Amat, nbc_Amat, &
                                  sca_ymat, desc_sca_ymat, &
                                  als_nu_max, als_block_partition, &
                                  merge(lambda_krr, als_ridge_k, als_precond_type==als_precond_flat), &
                                  als_rho, als_nsteps, als_tol, &
                                  als_alpha_method, als_precond_type, &
                                  als_nnls_alpha, als_nnls_mode, &
                                  context, nprow, npcol, myrow, mycol, &
                                  w_params(:,1), alpha_als, als_info)
        if (als_info /= 0 .and. myrow == 0 .and. mycol == 0) then
          write(6, '("ALS-Ridge ScaLAPACK: solver returned info = ", i8)') als_info
        end if
        ! Store the learned alpha_nu in the module variable for later use
        if (allocated(als_alpha)) deallocate(als_alpha)
        allocate(als_alpha(als_nu_max))
        als_alpha(:) = alpha_als(:)
        if (myrow == 0 .and. mycol == 0) then
          do inu = 1, als_nu_max
            write(6, '("ALS-Ridge ScaLAPACK: alpha(", i4, ") = ", ES15.6E3)') inu, als_alpha(inu)
          end do
        end if
        deallocate(alpha_als)
      end block

    case(fit_online_svd)
      ! Should not reach here — online fit has its own entry point in mld.F90
      call log_critical("fit_online_svd (mld_fit_type=6) must not use the standard ScaLAPACK parameter path")
      stop 'FATAL: fit_online_svd reached train_snap_get_parameters_linear_scalapack'

  end select 
  
  do ii = 1, size(w_params,1)
    call pdelset(sca_w_params, ii, 1, desc_sca_w_params, w_params(ii,1))
  end do 

  !$! ! ToScaRemove
  !$! if (rangml == 0 ) then
  !$!   do ii = 1, size(w_params,1)
  !$!     write(321, '(ES25.15E3)') w_params(ii,1)
  !$!   end do
  !$! end if
  !$! !end ToScaRemove
  call blacs_barrier(context, 'A')
  
  _MLD_END_
end subroutine train_snap_get_parameters_linear_scalapack


subroutine train_snap_get_parameters_kernel_scalapack()

  use ml_in_ndm_module, only: rangml, &
                              mld_fit_type, &
                              fit_home_hb, fit_lapack_qr, &
                              fit_lapack_ortho, fit_lapack_svd, fit_online_svd, svd_rcond

  use module_fit_ScaMatrix, only: prepare_sca_Amat_big, scalapack_lsystem_by_homeLU, & 
                                  scalapack_lsystem_by_QR, & 
                                  scalapack_square_lsystem_by_cholesky, &
                                  scalapack_lsystem_by_SVD

  use module_ml_scalapack, only: context, nbr_Amat, dimr_sca_Amat, &
                                 dimr_sca_Amat_big, dimc_sca_Amat_big, &
                                 sca_Amat_big, desc_sca_Amat_big, &
                                 nbr_Amat_big, nbc_Amat_big,  & 
                                 dimr_sca_Amat_big, dimc_sca_Amat_big, &
                                 sca_ymat_qr_svd, desc_sca_ymat_qr_svd, dimr_sca_ymat_qr_svd, &
                                 dimc_sca_ymat_qr_svd, &
                                 sca_w_params, desc_sca_w_params, dimr_sca_w_params, dimc_sca_w_params, &
                                 l_dimc_sca_w_params, l_dimr_sca_w_params, nbr_w_params, nbc_w_params,  &
                                sca_ymat, desc_sca_ymat, dimr_sca_ymat,& 
                                myrow, mycol, nprow, npcol, iam, context
  use module_scalapack_tools, only:   allocation_scalapack_matrix

  use snap, only: w_params, ymat
  use mld_logger

  implicit none
  integer :: rank_sca
  integer :: ii
  _NAMECURRENT_("train_snap_get_parameters_kernel_scalapack")



  _MLD_BEGIN_

  ! ToScaRemove
  if (rangml==0) write(6,*) 'Before prepare_sca_Amat_big'
  ! end ToScaRemove
  ! sca_w_params sca version for w_params
  !         |      ...   |
  !         |      ...   |
  !       = |  f(D) x 1  |  f(D) is line numbers of Amat
  !         |      ...   |
  !         |      ...   |
  !TODO WTF was nbr_Amat 
  nbr_w_params = nbr_Amat 
  nbc_w_params = 1
  !TODO WTF was dimr_sca_Amat 
  !if (rangml ==0) write(*,*) 'w_params dream', dimr_sca_Amat, size(w_params,1)
  dimr_sca_w_params = dimr_sca_Amat
  dimc_sca_w_params = 1
  call allocation_scalapack_matrix (" sca_w_params ", sca_w_params, desc_sca_w_params, &
                                  dimr_sca_w_params, dimc_sca_w_params, &
                                  l_dimr_sca_w_params, l_dimc_sca_w_params, nbr_w_params, nbc_w_params, & 
                                  myrow, mycol, nprow, npcol, context, iam  )

  ! set sca_ymat from ymat ... 
  do ii = 1, dimr_sca_ymat
    call pdelset(sca_ymat, ii, 1, desc_sca_ymat, ymat(ii,1)) 
    !debug write(24,'(2i8, es20.10)') iam, ii, ymat(ii,1)
  end do 
  call blacs_barrier(context, 'A')

  call prepare_sca_Amat_big
  call blacs_barrier(context, 'A')

  ! ToScaRemove
     if (rangml==0) write(6,*) 'Out from  prepare_sca_Amat_big'
  ! end ToScaRemove



  select case(mld_fit_type)

    !this is equivalent to HOME - only square matrix - Working tested. 
    case(fit_home_hb)

     call scalapack_lsystem_by_homeLU(dimr_sca_Amat_big, dimc_sca_Amat_big, sca_Amat_big, desc_sca_Amat_big, &
                                             nbr_Amat_big, nbc_Amat_big, &
                                             dimr_sca_ymat_qr_svd, dimc_sca_ymat_qr_svd, &
                                             sca_ymat_qr_svd, desc_sca_ymat_qr_svd, w_params)

    !this is QR -  general MxN matrix but full rank - Working tested.
    case (fit_lapack_qr)
      call scalapack_lsystem_by_QR(dimr_sca_Amat_big, dimc_sca_Amat_big, sca_Amat_big, desc_sca_Amat_big, & 
                                  nbr_Amat_big, nbc_Amat_big, &
                                  dimr_sca_ymat_qr_svd, dimc_sca_ymat_qr_svd, &
                                  sca_ymat_qr_svd, desc_sca_ymat_qr_svd, w_params)

    !this is Cholesky - only square matrix positive definite. 
    case (fit_lapack_ortho)
      call scalapack_square_lsystem_by_cholesky (dimr_sca_Amat_big, dimc_sca_Amat_big, sca_Amat_big, desc_sca_Amat_big, & 
                                  nbr_Amat_big, nbc_Amat_big, &
                                 dimr_sca_ymat_qr_svd, dimc_sca_ymat_qr_svd, &
                                 sca_ymat_qr_svd, desc_sca_ymat_qr_svd, w_params)
  
    !this is home made SVD for LSS problems, general matrix
    case(fit_lapack_svd)
      call scalapack_lsystem_by_SVD (dimr_sca_Amat_big, dimc_sca_Amat_big, sca_Amat_big, desc_sca_Amat_big, & 
                           nbr_Amat_big, nbc_Amat_big, &
                           dimr_sca_ymat_qr_svd, dimc_sca_ymat_qr_svd, &
                           sca_ymat_qr_svd, desc_sca_ymat_qr_svd, w_params, &
                           svd_rcond, rank_sca, context) 

    case(fit_online_svd)
      ! Should not reach here — online fit has its own entry point in mld.F90
      call log_critical("fit_online_svd (mld_fit_type=6) must not use the standard kernel ScaLAPACK path")
      stop 'FATAL: fit_online_svd reached train_snap_get_parameters_kernel_scalapack'

  end select 
  
  do ii = 1, size(w_params,1)
    call pdelset(sca_w_params, ii, 1, desc_sca_w_params, w_params(ii,1))
  end do 

  call blacs_barrier(context, 'A')
  
  _MLD_END_
end subroutine train_snap_get_parameters_kernel_scalapack



subroutine train_snap_get_parameters()

  use ml_in_ndm_module, only: mld_order,  &
                              mld_linear, mld_linear_extended, &
                              mld_quadratic, mld_type_quadratic, &
                              mld_type_quadratic_qml, mld_type_quadratic_qnml, mld_type_quadratic_bilinear, &
                              mld_polyc, mld_kernel, mld_type_quadratic_zaxa, & 
                              mld_type_quadratic_ZX

  use temporary_data_cov, only: dim_xdesc_patch
  use module_write_parameters, only: write_snap_parameters
  use module_ml_scalapack, only: scalapack_driver, context, nbr_Amat, dimr_sca_Amat, &
                                 sca_w_params, desc_sca_w_params, dimr_sca_w_params, dimc_sca_w_params, &
                                 l_dimc_sca_w_params, l_dimr_sca_w_params, nbr_w_params, nbc_w_params,  &
                                 myrow, mycol, nprow, npcol, iam, context
  use module_scalapack_tools, only:   allocation_scalapack_matrix

  use snap, only: w_params, dim_xdesc_linear, dim_xdesc_full
use mld_logger

  implicit none
  integer :: ii 
_NAMECURRENT_("train_snap_get_parameters")


  _MLD_BEGIN_

  if (mld_order == mld_linear) then
    !if ((.not.weighted).and.(.not.weighted_3ch)) then
    if (scalapack_driver) then
      call train_snap_get_parameters_linear_scalapack()
    else
    call train_snap_get_parameters_linear
    end if 
    !call train_snap_get_parameters_kernel
    !end if
    !if (weighted.and.(.not.(weighted_3ch))) then
    ! call  train_snap_get_parameters_precondition_inversion_matrix(dim_xdesc/2, dim_xdesc+1)
    !end if
    !if (weighted.and.weighted_3ch) then
    ! call  train_snap_get_parameters_precondition_inversion_matrix(dim_xdesc/3, dim_xdesc+1)
    !end if
  end if

  if (mld_order == mld_linear_extended) then
    call init_objective_minimization
    call train_snap_get_parameters_quadratic_objective_minimization
  end if


  if (mld_order == mld_quadratic) then

    ! Full quadratic version on E ...
    if (mld_type_quadratic == mld_type_quadratic_qml) then
      if (scalapack_driver) then 
        call train_snap_get_parameters_linear_scalapack()
      else
        call train_snap_get_parameters_linear
      end if
    end if

    ! The quadratic version on E - E_LML
    if ((mld_type_quadratic == mld_type_quadratic_qnml).or. & 
        (mld_type_quadratic == mld_type_quadratic_zaxa).or. &
        (mld_type_quadratic == mld_type_quadratic_ZX)) then
      !call train_snap_get_parameters_precondition_inversion_matrix(dim_xdesc, dim_xdesc_quadratic)
      call train_snap_get_parameters_precondition_inversion_matrix(dim_xdesc_linear, dim_xdesc_patch, dim_xdesc_full)
    end if

    if (mld_type_quadratic == mld_type_quadratic_bilinear) then
      call train_snap_get_parameters_quadratic_bilinear_inversion_matrix
    end if

  end if

  if (mld_order == mld_polyc) then
    ! The polyc version on E - E_LML
    call train_snap_get_parameters_precondition_inversion_matrix(dim_xdesc_linear, dim_xdesc_patch, dim_xdesc_full)
  end if

  if (mld_order == mld_kernel) then
    ! The kernel version on E - E_LML
    call train_snap_get_parameters_precondition_inversion_matrix(dim_xdesc_linear, dim_xdesc_patch, dim_xdesc_full)

    !if (scalapack_driver) then 
    !  call train_snap_get_parameters_linear_scalapack()
    !else
    !call train_snap_get_parameters_linear
    !end if 


  end if
  
  if (scalapack_driver)  then 
  
    call blacs_barrier(context, 'A')
    if (allocated(sca_w_params)) deallocate(sca_w_params)
    if (allocated(desc_sca_w_params)) deallocate(desc_sca_w_params)
    nbr_w_params = nbr_Amat 
    nbc_w_params = 1
    dimr_sca_w_params = dimr_sca_Amat
    dimc_sca_w_params = 1
    call allocation_scalapack_matrix (" sca_w_params ", sca_w_params, desc_sca_w_params, &
                                  dimr_sca_w_params, dimc_sca_w_params, &
                                  l_dimr_sca_w_params, l_dimc_sca_w_params, nbr_w_params, nbc_w_params, & 
                                  myrow, mycol, nprow, npcol, context, iam  )
    do ii = 1, size(w_params,1)
      call pdelset(sca_w_params, ii, 1, desc_sca_w_params, w_params(ii,1))
    end do 
                                
  end if 

  !if (iam==0) write(6,*) 'before write_snap_parameters' 
  ! TOLD call write_snap_parameters
  if (scalapack_driver)  call blacs_barrier(context, 'A')
  !if (iam==0) write(6,*) 'after write_snap_parameters' 


  _MLD_END_
end subroutine train_snap_get_parameters


subroutine split_Amat_in_Lmat_Qmat(dim_descriptor_base, dim_patch, dim_descriptor_full)

  !     --------- M columns ----------
  !    | d_1(1)             d_1(M)    |     |
  !    |   .                  .       |     |
  ! A =|   .        ...       .       |    f(D) lines
  !    |   .                  .       |     |
  !    | d_f(D)(1)          d_f(D)(M) |     |
  !
  !    will be decomposed into two matrix: LMAT and QMAT
  !
  !            --------- M columns ----------
  !         |                                |
  !  LMAT = |                                |  D + Dpatch  + 1 lines
  !         |                                |
  !
  !
  !            --------- M columns ----------
  !         |                                |
  !  QMAT = |                                |  f(D) - D - Dpatch - 1  lines
  !         |                                |
  !  the only assumption is the knowledge of
  !  dim_desciptor_base = D,  dim_descriptor_full=f(D), Dpatch= dimension of the patch 


  use ml_in_ndm_module, only: rangml
  !use temporary_data_cov, only: dim_xdesc
  use snap, only: Amat
  use module_mld_quadratic, only: Lmat, Qmat
  use module_ml_scalapack, only: dimr_sca_Amat, dimc_sca_Amat, sca_Amat, desc_sca_Amat, nbr_Amat, nbc_Amat, &
                                 sca_Lmat, desc_sca_Lmat, &
                                 dimr_sca_Lmat, dimc_sca_Lmat, &
                                 l_dimr_sca_Lmat, l_dimc_sca_Lmat, nbr_Lmat, nbc_Lmat, & 
                                 sca_Qmat, desc_sca_Qmat, &
                                 dimr_sca_Qmat, dimc_sca_Qmat, &
                                 l_dimr_sca_Qmat, l_dimc_sca_Qmat, nbr_Qmat, nbc_Qmat, & 
                                 myrow, mycol, nprow, npcol, context, iam, & 
                                 scalapack_driver, context
                                 
  use module_scalapack_tools, only: allocation_scalapack_matrix

  use mld_logger

  implicit none

  integer, intent(in)  :: dim_descriptor_base,  dim_patch, dim_descriptor_full
  integer  :: size_Aline

  _NAMECURRENT_("split_Amat_in_Lmat_Qmat")



  _MLD_BEGIN_
  ! The version E - E_LML
  if (scalapack_driver) then
    size_Aline = dimr_sca_Amat 
  else 
    size_Aline = size(Amat, 1)
  end if 
  if (size_Aline /= dim_descriptor_full) then
    if (rangml == 0) write (6, *) 'size_Aline vs dim_descriptor', size_Aline, dim_descriptor_full
    stop 'huge problem in split_Amat_in_Lmat_Qmat size_Aline vs dim_descriptor_full'
  end if

  if (scalapack_driver) then
    ! allocate desc for sca_Lmat
    nbr_Lmat = nbr_Amat 
    nbc_Lmat = nbc_Amat
    dimr_sca_Lmat = dim_descriptor_base + dim_patch 
    dimc_sca_Lmat = dimc_sca_Amat
    call allocation_scalapack_matrix (" sca_Lmat ", sca_Lmat, desc_sca_Lmat, &
                                      dimr_sca_Lmat, dimc_sca_Lmat, &
                                      l_dimr_sca_Lmat, l_dimc_sca_Lmat, nbr_Lmat, nbc_Lmat, & 
                                      myrow, mycol, nprow, npcol, context, iam  )
    call blacs_barrier(context, 'A')
    nbr_Qmat = nbr_Amat 
    nbc_Qmat = nbc_Amat
    dimr_sca_Qmat = size_Aline - dim_descriptor_base - dim_patch 
    dimc_sca_Qmat = dimc_sca_Amat
    call allocation_scalapack_matrix (" sca_Qmat ", sca_Qmat, desc_sca_Qmat, &
                                      dimr_sca_Qmat, dimc_sca_Qmat, &
                                      l_dimr_sca_Qmat, l_dimc_sca_Qmat, nbr_Qmat, nbc_Qmat, & 
                                      myrow, mycol, nprow, npcol, context, iam  )
    call blacs_barrier(context, 'A')
    
    ! copy a into b: in this care the first part of sca_Amat into sca_Lmat
    !call pdgemr2d(m, n, a, ia, ja, desca, b, ib, jb, descb, ictxt)
    call pdgemr2d(dimr_sca_Lmat, dimc_sca_Lmat, sca_Amat, 1, 1, desc_sca_Amat, &
                           sca_Lmat, 1, 1, desc_sca_Lmat, context)

    ! copy a into b: in this care the second part of sca_Amat into sca_Qmat 
    !call pdgemr2d(m, n, a, ia, ja, desca, b, ib, jb, descb, ictxt)
    call pdgemr2d(dimr_sca_Qmat, dimc_sca_Qmat, sca_Amat, dimr_sca_Lmat + 1, 1, desc_sca_Amat, &
                            sca_Qmat, 1, 1, desc_sca_Qmat, context)
  else 
    ! size(Amat,1) -> f(D)
    ! size(Amat,2) -> M
  
    if (size_Aline <= (dim_descriptor_base + dim_patch) ) then
      if (rangml == 0) write (6, *) 'size_Aline lower than D ', size_Aline, dim_descriptor_full
      stop 'huge problem in split_Amat_in_Lmat_Qmat size_Aline vs dim_descriptor_base'
    end if
  
    if (allocated(Lmat)) deallocate (Lmat); allocate (Lmat(dim_descriptor_base + dim_patch, size(Amat, 2)))
    if (allocated(Qmat)) deallocate (Qmat); allocate (Qmat(size_Aline - dim_descriptor_base - dim_patch, size(Amat, 2)))
    Lmat(1:dim_descriptor_base + dim_patch, :) = Amat(1:dim_descriptor_base + dim_patch, :)
    Qmat(1:size_Aline - dim_descriptor_base - dim_patch, :) = Amat(dim_descriptor_base + dim_patch + 1:size_Aline, :)
  end if 

  _MLD_END_

end subroutine split_Amat_in_Lmat_Qmat



subroutine precond_linear_from_Lmat
  ! the purposes of this subroputine: 
  !    - get the  linear parameters from Lmat using train_snap_get_parameters_linear
  !    - make a copy of (sca_)ymat into (sca_)ymat_copy
  !    - redefine (sca_)ymat by  (sca_)ymat - LML prediction
  ! is patched using scalapack_driver 
  use ml_in_ndm_module, only: rangml
  use snap, only: Amat, w_params, ymat, fit_snap
  use module_mld_quadratic, only: Lmat
  use derived_types, only: config_real
  use module_evaluate_parameters, only: product_w_params_Amat
  use module_ml_scalapack, only: dimr_sca_Amat, dimc_sca_Amat, sca_Amat, desc_sca_Amat, &
                                 nbr_Amat, nbc_Amat, l_dimr_sca_Amat, l_dimc_sca_Amat, &
                                 sca_Lmat, desc_sca_Lmat, &
                                 sca_ymat, desc_sca_ymat, & 
                                 dimr_sca_Lmat, dimc_sca_Lmat, &
                                 nbr_Lmat, nbc_Lmat, & 
                                 myrow, mycol, nprow, npcol, context, iam, & 
                                 scalapack_driver
  use module_scalapack_tools, only: allocation_scalapack_matrix
  use mld_logger

  implicit none

  integer  :: it, Ncolm, nounit1, nounit2, nounit3
  real(kind=kind(0.d0))      :: term

  _NAMECURRENT_("precond_linear_from_Lmat")

  _MLD_BEGIN_
  if (scalapack_driver) then
    if (allocated(w_params)) deallocate (w_params); allocate (w_params(dimr_sca_Lmat, 1))

    !reallocate sca_Amat with sca_Lmat shape.
    if (allocated(sca_Amat)) deallocate (sca_Amat)
    if (allocated(desc_sca_Amat)) deallocate (desc_sca_Amat)
    nbr_Amat = nbr_Lmat 
    nbc_Amat = nbc_Lmat
    dimr_sca_Amat = dimr_sca_Lmat 
    dimc_sca_Amat = dimc_sca_Lmat
    call allocation_scalapack_matrix (" sca_Amat_from_Lmat ", sca_Amat, desc_sca_Amat, &
                                      dimr_sca_Amat, dimc_sca_Amat, &
                                      l_dimr_sca_Amat, l_dimc_sca_Amat, nbr_Amat, nbc_Amat, & 
                                      myrow, mycol, nprow, npcol, context, iam  )

    ! copy sca_Lmat into sca_Amat
    call pdgemr2d(dimr_sca_Lmat, dimc_sca_Lmat, sca_Lmat, 1, 1, desc_sca_Lmat, &
                                      sca_Amat, 1, 1, desc_sca_Amat, context)

    call train_snap_get_parameters_linear_scalapack
    Ncolm = dimc_sca_Amat
  else 
    if (allocated(Amat)) deallocate (Amat); allocate (Amat(size(Lmat, 1), size(Lmat, 2)))
    Amat(:, :) = Lmat(:, :)
    if (allocated(w_params)) deallocate (w_params); allocate (w_params(size(Lmat, 1), 1))
    call train_snap_get_parameters_linear
    Ncolm = size(Amat,2)
  end if !scalapack_driver 

  if (rangml==0) open (file="energy_deviation_linear.milady", newunit=nounit1, action='write')
  if (rangml==0) open (file="force_deviation_linear.milady", newunit=nounit2, action='write')
  if (rangml==0) open (file="stress_deviation_linear.milady", newunit=nounit3, action='write')
  do it = 1, Ncolm 
    !term = DOT_PRODUCT(w_params(:, 1), Amat(:, it))
    call product_w_params_Amat(term, it)
    if (rangml==0) then 
      if (fit_snap(it)%energy) write (nounit1, '(i9," ",3e20.10, " ", (a), " ", (a))') it, ymat(it, 1) - term, term, ymat(it, 1), config_real(fit_snap(it)%iconf)%class, trim(config_real(fit_snap(it)%iconf)%filename)
      if (fit_snap(it)%force) write (nounit2, '(i9," ",3e20.10, " ", (a), " ", (a))') it, ymat(it, 1) - term, term, ymat(it, 1), config_real(fit_snap(it)%iconf)%class, trim(config_real(fit_snap(it)%iconf)%filename)
      if (fit_snap(it)%stress) write (nounit3, '(i9," ",3e20.10, " ", (a), " ", (a))') it, ymat(it, 1) - term, term, ymat(it, 1), config_real(fit_snap(it)%iconf)%class, trim(config_real(fit_snap(it)%iconf)%filename)
    end if

    ymat(it, 1) = ymat(it, 1) - term
    if (scalapack_driver) then
      call pdelset(sca_ymat, it, 1, desc_sca_ymat,  ymat(it, 1))
    end if 
  end do

  if (rangml==0) close (nounit1, status='keep')
  if (rangml==0) close (nounit2, status='keep')
  if (rangml==0) close (nounit3, status='keep')

  _MLD_END_

end subroutine precond_linear_from_Lmat


subroutine precond_nonlinear_from_Qmat
  !   - reallocate w_params 
  !   - copy (sca_)Qmat into (sca_)Amat 
  !   - reallocate w_paramsfor Qmat part. (sca_w_params) are reallocated elsewhere. 
  !   - get the parameters for Qmat part
  use module_kind_variables, only: kind_double
  use ml_in_ndm_module, only: rangml, mld_polyc,  &
                              mld_order, mld_polyc, mld_quadratic, &
                              mld_kernel, mld_linear

  use derived_types, only: config_real
  use snap, only: Amat, w_params, ymat, fit_snap
  use module_mld_quadratic, only: Qmat
  use module_evaluate_parameters, only: product_w_params_Amat
  use module_ml_scalapack, only: sca_Amat, desc_sca_Amat, nbr_Amat, nbc_Amat, &
                                 dimr_sca_Amat, dimc_sca_Amat, &
                                 l_dimr_sca_Amat, l_dimc_sca_Amat, &
                                 nbr_Qmat, nbc_Qmat, dimr_sca_Qmat, dimc_sca_Qmat, &
                                 sca_Qmat, desc_sca_Qmat, &
                                  myrow, mycol, nprow, npcol, context, iam, & 
                                  scalapack_driver
                                  
                                  
                                  

  use module_scalapack_tools, only: allocation_scalapack_matrix
  use module_variance, only: sigmam_e, sigmam_f, sigmam_s

  implicit none 
  real(kind_double)     :: term, delta_corr
  integer :: it, rnumber,  nounit92, nounit93, nounit94
  integer :: count_e, count_f, count_s

  ! --------------------- Amat here is Qmat ------------------------!

  if (scalapack_driver) then
    if (allocated(w_params)) deallocate (w_params); allocate (w_params(dimr_sca_Qmat, 1))
    if (allocated(sca_Amat)) deallocate (sca_Amat)
    if (allocated(desc_sca_Amat)) deallocate (desc_sca_Amat)
    nbr_Amat = nbr_Qmat 
    nbc_Amat = nbc_Qmat
    dimr_sca_Amat = dimr_sca_Qmat 
    dimc_sca_Amat = dimc_sca_Qmat
    call allocation_scalapack_matrix (" sca_Amat_from_Qmat ", sca_Amat, desc_sca_Amat, &
                                        dimr_sca_Amat, dimc_sca_Amat, &
                                        l_dimr_sca_Amat, l_dimc_sca_Amat, nbr_Amat, nbc_Amat, & 
                                        myrow, mycol, nprow, npcol, context, iam  )
    call blacs_barrier(context, 'A')
    call pdgemr2d(dimr_sca_Qmat, dimc_sca_Qmat, sca_Qmat, 1, 1, desc_sca_Qmat, &
                                        sca_Amat, 1, 1, desc_sca_Amat, context)

    if (mld_order == mld_kernel) then 
      call train_snap_get_parameters_kernel_scalapack
      !call train_snap_get_parameters_linear_scalapack
    else
      call train_snap_get_parameters_linear_scalapack
    end if     
  else 
    if (allocated(w_params)) deallocate (w_params); allocate (w_params(size(Qmat, 1), 1))
    deallocate (Amat); allocate (Amat(size(Qmat, 1), size(Qmat, 2)))
    Amat(:, :) = Qmat(:, :)
    if (mld_order == mld_kernel) then
      !call train_snap_get_parameters_linear
      call train_snap_get_parameters_kernel
    else
      call train_snap_get_parameters_linear
    end if
    
  end if 
  !-----------------------------------------------------------------
  if (rangml==0) then 
    select case (mld_order)
      case (mld_polyc)
        open (file="energy_deviation_polyc.milady", newunit=nounit92, action='write')
        open (file="force_deviation_polyc.milady", newunit=nounit93, action='write')
        open (file="stress_deviation_polyc.milady", newunit=nounit94, action='write')
  
      case (mld_quadratic)
        open (file="energy_deviation_quadratic.milady", newunit=nounit92, action='write')
        open (file="force_deviation_quadratic.milady", newunit=nounit93, action='write')
        open (file="stress_deviation_quadratic.milady", newunit=nounit94, action='write')
  
  
      case (mld_kernel)
        open (file="energy_deviation_kernel.milady", newunit=nounit92, action='write')
        open (file="force_deviation_kernel.milady", newunit=nounit93, action='write')
        open (file="stress_deviation_kernel.milady", newunit=nounit94, action='write')
  
      case (mld_linear)
        open (file="energy_deviation_ch.milady", newunit=nounit92, action='write')
        open (file="force_deviation_ch.milady", newunit=nounit93, action='write')
        open (file="stress_deviation_ch.milady", newunit=nounit94, action='write')
  
      case default
        if (rangml == 0) write (6, *) 'No precondition by matrix for this potential type, mld_order', mld_order
        stop 'wrorng precondition for this potential type'
    end select
  end if

  if (scalapack_driver) then
    rnumber = dimc_sca_Amat 
  else 
    rnumber = size(Amat, 2)
  end if

  sigmam_e = 0.d0
  sigmam_f = 0.d0
  sigmam_s = 0.d0
  count_e = 0 
  count_f = 0 
  count_s = 0 

  do it = 1, rnumber
    !term = DOT_PRODUCT(w_params(:, 1), Amat(:, it))
    call product_w_params_Amat(term, it)
    !term = 0 
    delta_corr = ymat(it, 1) - term

    !if (rangml==0) then   
      if (fit_snap(it)%energy) then
    if (rangml==0) then   
          write (nounit92, '(i9," ",3e20.10, " ", (a), " ", (a))') it, delta_corr, term, ymat(it, 1), config_real(fit_snap(it)%iconf)%class, trim(config_real(fit_snap(it)%iconf)%filename)
    end if 
        sigmam_e = sigmam_e + delta_corr**2
        count_e = count_e+1
      end if
      if (fit_snap(it)%force) then 
        if (rangml==0) then
          write (nounit93, '(i9," ",3e20.10, " ", (a), " ", (a))') it, delta_corr, term, ymat(it, 1), config_real(fit_snap(it)%iconf)%class, trim(config_real(fit_snap(it)%iconf)%filename)
        end if
        sigmam_f = sigmam_f + delta_corr**2
        count_f = count_f+1
      end if
      if (fit_snap(it)%stress) then
        if (rangml==0) then
          write (nounit94, '(i9," ",3e20.10, " ", (a), " ", (a))') it, delta_corr, term, ymat(it, 1), config_real(fit_snap(it)%iconf)%class, trim(config_real(fit_snap(it)%iconf)%filename)
        end if
        sigmam_s = sigmam_s + delta_corr**2
        count_s = count_s+1
      end if    
    !end if 
  end do
  
  sigmam_e=sigmam_e/dble(count_e-1)
  sigmam_f=sigmam_f/dble(count_f-1)
  sigmam_s=sigmam_s/dble(count_s-1)
  !write (6, *) '!!! rangml, sigmam_f !!!', rangml, sigmam_f

  if (rangml==0) close (nounit92, status='keep')
  if (rangml==0) close (nounit93, status='keep')
  if (rangml==0) close (nounit94, status='keep')
 
end  subroutine precond_nonlinear_from_Qmat 


subroutine train_snap_get_parameters_precondition_inversion_matrix(dim_descriptor_base, dim_patch, dim_descriptor_full)
  ! whar this subroutine do (serial and parallel version)
  !  1- split Amat = Lmat + Qmat by calling split_Amat_in_Lmat_and_Qmat
  !   - make a copy of (sca_ymat) into (sca_)ymat_copy
  !   - allocate  global w_params_final
  !  2- copy linear part into w_params_final  
  !   - copy (sca_)Lmat into (sca_)Amat 
  !   - get the paramteters for linear part 
  !  3- get the parameters for Qmat part by a call precond_nonlinear_from_Qmat
  !  4- compute the deviantion LML from Qmat
  !  5- final touch: reput everything in place
  !   - complete w_params_final with Qmat part 
  !   - reallocate w_params and copy w_params_final into w_params 
  !   - fill sca_w_params with w_params 
  !   - reallocate initial (sca_)Amat and fill that with (sca_)Lmat and (sca_)Qmat
  !   - deallocate (sca_)Lmat and (sca_) Qmat   

  use ml_in_ndm_module, only: rangml

  use snap, only: Amat, w_params, ymat
  use module_mld_quadratic, only: Lmat, Qmat, ymat_copy, w_params_final
  use module_evaluate_parameters, only: product_w_params_Amat
  use module_ml_scalapack, only: scalapack_driver, dimr_sca_Amat, &
                                 sca_ymat_copy, desc_sca_ymat_copy, &
!                                 dimr_sca_ymat_copy,  &
!                                  l_dimc_sca_ymat_copy, &
!                                 nbr_ymat_copy,  & 
!                                 nbr_ymat, & 
!                                 dimr_sca_ymat,  &
                                 sca_ymat, desc_sca_ymat, &
                                 nbr_Amat, nbc_Amat, & 
                                 dimr_sca_Amat, dimc_sca_Amat, &
                                 l_dimr_sca_Amat, l_dimc_sca_Amat, &
                                 sca_Amat, desc_sca_Amat, &
!                                 nbr_Qmat,  & 
                                 dimr_sca_Qmat, dimc_sca_Qmat, &
!                                 l_dimr_sca_Qmat,  &
                                 sca_Qmat, desc_sca_Qmat, &
                                 nbr_Lmat, nbc_Lmat, dimr_sca_Lmat, dimc_sca_Lmat, &
                                 sca_Lmat, desc_sca_Lmat , & 
                                 sca_w_params, desc_sca_w_params, dimr_sca_w_params, &
                                 dimc_sca_w_params, l_dimr_sca_w_params, l_dimc_sca_w_params, &
                                 nbr_w_params, nbc_w_params, &
                                 myrow, mycol, nprow, npcol, context, iam
  use module_scalapack_tools, only: allocation_scalapack_matrix   
  use mld_logger

  implicit none

  integer, intent(in)  :: dim_descriptor_full, dim_patch, dim_descriptor_base
  integer  :: size_Aline, ii
!  real(kind(1.d0))     :: term

  _NAMECURRENT_("train_snap_get_parameters_precondition_inversion_matrix")


  _MLD_BEGIN_
  ! The version E - E_LML
  ! size(Amat,1) -> f(D) -> dim_descriptor_full
  ! size(Amat,2) -> M

  if (scalapack_driver) then
    size_Aline = dimr_sca_Amat
  else 
    size_Aline = size(Amat, 1)
  end if
  
  ! dim_descriptor_base = dim_xdesc
  if (dim_descriptor_full /= size_Aline ) then
    if (rangml == 0) write (6, *) 'ML: problems in dimensions in '//NAMECURRENT , &
    dim_descriptor_full, size_Aline
    stop
  end if
  call split_Amat_in_Lmat_Qmat(dim_descriptor_base, dim_patch,  dim_descriptor_full)

  !Ymat should be copied ...
  if (scalapack_driver) then
    !$$$! nbr_ymat_copy = nbr_ymat
    !$$$! nbc_ymat_copy = nbc_ymat
    !$$$! dimr_sca_ymat_copy = dimr_sca_ymat
    !$$$! dimc_sca_ymat_copy = dimc_sca_ymat
    !$$$! call allocation_scalapack_matrix (" sca_ymat_copy ", sca_ymat_copy, desc_sca_ymat_copy, &
    !$$$!                                   dimr_sca_ymat_copy, dimc_sca_ymat_copy, &
    !$$$!                                   l_dimr_sca_ymat_copy, l_dimc_sca_ymat_copy, nbr_ymat_copy, nbc_ymat_copy, & 
    !$$$!                                   myrow, mycol, nprow, npcol, context, iam  )
    !$$$! ! sca_ymat is copied into sca_ymat_copy
    !$$$! call pdgemr2d(dimr_sca_ymat, dimc_sca_ymat, sca_ymat, 1, 1, desc_sca_ymat, &
    !$$$!                                   sca_ymat_copy, 1, 1, desc_sca_ymat_copy, context)

    ! ToScaRemove? - maybe to have a global copy is not a bad ideea.
    if (allocated(ymat_copy)) deallocate (ymat_copy); allocate(ymat_copy(dimc_sca_Amat, 1))
    ymat_copy = ymat

  else 
    if (allocated(ymat_copy)) deallocate (ymat_copy); allocate (ymat_copy(size(Amat, 2), 1))
    ymat_copy = ymat
  end if 

  ! this subroutine get the parameters and the new ymat with the diff E_DFT - E_LML
  call precond_linear_from_Lmat

  if (allocated(w_params_final)) deallocate (w_params_final); allocate (w_params_final(dim_descriptor_full, 1))
  w_params_final(1:dim_descriptor_base + dim_patch, 1) = w_params(1:dim_descriptor_base + dim_patch, 1)

  call precond_nonlinear_from_Qmat


  ! final touch ...
  if (scalapack_driver) then
    size_Aline = dimr_sca_Qmat
  else
    size_Aline = size(Amat,1) 
  end if 

  ! recopy w_params ...
  w_params_final(dim_descriptor_base + dim_patch + 1:dim_descriptor_full, 1) = w_params(1:size_Aline,1)
  if (allocated(w_params)) deallocate (w_params); allocate (w_params(dim_descriptor_full, 1))
  w_params(1:dim_descriptor_full, :) = w_params_final(1:dim_descriptor_full, :)
  deallocate (w_params_final)

  
  if (scalapack_driver) then 

   !$$$!call pdgemr2d(dimr_sca_ymat, dimc_sca_ymat, sca_ymat_copy, 1, 1, desc_sca_ymat_copy, &
   !$$$!                                   sca_ymat, 1, 1, desc_sca_ymat, context)
    
    do ii = 1, size(ymat_copy,1)
      call pdelset(sca_ymat, ii, 1, desc_sca_ymat, ymat_copy(ii,1))
    end do 

    if (allocated(sca_ymat_copy))  deallocate(sca_ymat_copy)
    if (allocated(desc_sca_ymat_copy)) deallocate(desc_sca_ymat_copy)

    if (allocated(sca_Amat)) deallocate(sca_Amat)
    if (allocated(desc_sca_Amat)) deallocate(desc_sca_Amat)
    nbr_Amat = nbr_Lmat 
    nbc_Amat = nbc_Lmat 
    dimr_sca_Amat = dimr_sca_Qmat + dimr_sca_Lmat
    dimc_sca_Amat = dimc_sca_Qmat
    call allocation_scalapack_matrix (" sca_Amat_final ", sca_Amat, desc_sca_Amat, &
                                      dimr_sca_Amat, dimc_sca_Amat, &
                                      l_dimr_sca_Amat, l_dimc_sca_Amat, nbr_Amat, nbc_Amat, & 
                                      myrow, mycol, nprow, npcol, context, iam  )

    
    call pdgemr2d(dimr_sca_Lmat, dimc_sca_Lmat,  sca_Lmat, 1, 1, desc_sca_Lmat, &
                                      sca_Amat, 1, 1, desc_sca_Amat, context)

    call pdgemr2d(dimr_sca_Qmat, dimc_sca_Qmat,  sca_Qmat, 1, 1, desc_sca_Qmat, &
                                      sca_Amat, dimr_sca_Lmat+1, 1, desc_sca_Amat, context)
    if (allocated(sca_Qmat))   deallocate(sca_Qmat)
    if (allocated(desc_sca_Qmat)) deallocate(desc_sca_Qmat)
    if (allocated(sca_Lmat))   deallocate(sca_Lmat)
    if (allocated(desc_sca_Lmat)) deallocate(desc_sca_Lmat)


    if (allocated(sca_w_params)) deallocate(sca_w_params)
    if (allocated(desc_sca_w_params)) deallocate(desc_sca_w_params)
    nbr_w_params = nbr_Amat 
    nbc_w_params = 1
    dimr_sca_w_params = dimr_sca_Amat
    dimc_sca_w_params = 1
    call allocation_scalapack_matrix (" sca_w_params ", sca_w_params, desc_sca_w_params, &
                                  dimr_sca_w_params, dimc_sca_w_params, &
                                  l_dimr_sca_w_params, l_dimc_sca_w_params, nbr_w_params, nbc_w_params, & 
                                  myrow, mycol, nprow, npcol, context, iam  )

    do ii = 1, size(w_params,1)
      call pdelset(sca_w_params, ii, 1, desc_sca_w_params, w_params(ii,1))
    end do 

  else 
    ymat = ymat_copy
    deallocate(ymat_copy)
    deallocate (Amat); allocate (Amat(dim_descriptor_full, size(Qmat, 2)))
    Amat(1:dim_descriptor_base + dim_patch, :) = Lmat(1:dim_descriptor_base + dim_patch, :)
    deallocate (Lmat)
    Amat(dim_descriptor_base + dim_patch +1 : dim_descriptor_full, :) = Qmat(:, :)
    deallocate (Qmat)
  end if 

  _MLD_END_

end subroutine train_snap_get_parameters_precondition_inversion_matrix


subroutine train_snap_get_parameters_linear()

  use ml_in_ndm_module, only: rangml, svd_rcond, mld_fit_type, &
                              fit_home_hb, fit_lapack_qr, fit_lapack_qr_constraints, &
                              fit_lapack_ortho, fit_lapack_svd, fit_als, fit_online_svd, lambda_krr, &
                              descriptor_type, descriptor_body, descriptor_ace
  use module_input_als_fit, only: als_nsteps, als_tol, als_ridge_k, als_rho, &
                              als_alpha_method, als_precond_type, als_precond_flat, &
                              als_nnls_alpha, als_nnls_mode, &
                              als_nu_max, als_block_partition, als_alpha, &
                              set_als_block_partition
  use temporary_data_cov, only: dim_data_train
  use main_mld_mod       ! TODOcvw
  use snap, only: w_params, Amat, ymat, Bmat, zmat, i_fit_snap, weights_snap
  use mld_logger
  use math, only: dreal_matmul, dreal_inverse
  use module_serial_linear_solver, only: serial_lsystem_by_svd, serial_lsystem_by_ortho, &
                                         serial_lsystem_by_qr, serial_lsystem_by_constraints, &
                                         serial_lsystem_by_home
  use module_serial_als, only: serial_als_ridge_solve

  implicit none

  integer  :: i

  real(kind=kind(1.d0)), allocatable, dimension(:, :)      :: phi, Cmat, Cmat_transpose
  real(kind=kind(1.d0)), dimension(:, :), allocatable      :: ymat_qr_svd
  real(kind=kind(1.d0)) :: svd_rcond_local

  integer  ::  Mline, Ncolm

  _NAMECURRENT_("train_snap_get_parameters_linear")


  _MLD_BEGIN_
  !  Amat is the transpose of the design matrix, i.e. Amat =  P_descriptor(D) x  M_observations
  !                      M_obsv
  !           -------------------------------
  !         |                                |
  !   Amat  |                                |   is short-fat
  !         |                                |
  !           -------------------------------
  ! Amat is the transpose of the design matrix, i.e. Amat =  P_descriptor(D) x  M_observations X
  ! size(Amat,1) = dim_xdesc+1, dim_xdesc**2 + dim_xdesc +1, dim_kernel etc
  ! size(Amat,2) = number of observations.

  if (allocated(phi)) deallocate (phi); allocate (phi(size(Amat, 1), size(Amat, 1)))
  do i = 1, size(Amat,1) 
    phi(i,i) = 0.d0 
  end do 
  if (allocated(Cmat)) deallocate (Cmat); allocate (Cmat(size(Amat, 1), size(Amat, 2)))
  if (allocated(Cmat_transpose)) deallocate (Cmat_transpose); allocate (Cmat_transpose(size(Amat, 2), size(Amat, 1)))

  if (i_fit_snap /= dim_data_train) then
    if (rangml == 0) write (6, *) 'The problem is not correct sized', i_fit_snap, dim_data_train
    stop 'incorrect sized problem train_snap_get_parameters'
  end if

  Cmat(:, :) = Amat(:, :)
  ! the solution is (Amat W Amat^T)^-1 A W ymat or (Amat Cmat^T)^-1 Cmat ymat
  ! Cmat = Amat x W and W^T = W because diagonal.
  do i = 1, size(Cmat, 2)
    Cmat(:, i) = Amat(:, i)*weights_snap(i)
  end do

  Cmat_transpose(:, :) = transpose(Cmat)
  ! phi = matmul ( Amat, Cmat_transpose)
  call dreal_matmul(Amat, size(Amat, 1), size(Amat, 2), Cmat_transpose, size(Cmat_transpose, 1), size(Cmat_transpose, 2), &
                    phi, size(phi, 1), size(phi, 2))

  ! solving with regularization
  if (lambda_krr < 0) lambda_krr = 0.d0
  do i = 1, size(Amat, 1)
    phi(i,i) = phi(i,i) + lambda_krr
  end do

  Mline = size(phi, 1)
  Ncolm = size(phi, 2)
  ! prepare the ymat_qr_svd = Amat W ymat  = Cmat ymat
  if (allocated(ymat_qr_svd)) deallocate (ymat_qr_svd); allocate (ymat_qr_svd(size(Cmat, 1), size(ymat, 2)))
  call dreal_matmul(Cmat, size(Cmat, 1), size(Cmat, 2), ymat, size(ymat, 1), size(ymat, 2), &
                    ymat_qr_svd, size(ymat_qr_svd, 1), size(ymat_qr_svd, 2))


  !call prepare_serial_phi (Amat, weights_snap, lambda_krr, Cmat, phi, ymat_qr_svd)
  !Mline = size(phi, 1)
  !Ncolm = size(phi, 2)

  select case (mld_fit_type)


  case (fit_home_hb)

    !call serial_lsystem_by_home (phi, Cmat, ymat, w_params)
    call serial_lsystem_by_home (phi, ymat_qr_svd, Mline, Ncolm, w_params)

  case (fit_lapack_qr)
    ! QR based on Lapack dgels LLS general A MxN matrix.
    ! It is assumed that A has full rank.
    call serial_lsystem_by_qr(phi, ymat_qr_svd, Mline, Ncolm, w_params)

    if (allocated (ymat_qr_svd)) deallocate(ymat_qr_svd)
    if (allocated(phi)) deallocate(phi)

    !mld_fit_type=3 ortho decomposition, rank estimation ...
  case (fit_lapack_ortho)
    ! ortho solution based on dgelsy
    ! DGELSY computes the minimum-norm solution to a real linear least
    ! squares problem:
    ! minimize || A * X - B ||
    ! using a complete orthogonal factorization of A.  A is an M-by-N
    ! matrix which may be rank-deficient.
    call serial_lsystem_by_ortho (phi, ymat_qr_svd, Mline, Ncolm, w_params)
  
    if (allocated(ymat_qr_svd)) deallocate(ymat_qr_svd)
    if (allocated(phi)) deallocate(phi)

    !mld_fit_type=4 SVD, rank estimation
  case (fit_lapack_svd)
    ! svd solution based on dgelsd
    !DGELSD computes the minimum-norm solution to a real linear least
    !squares problem: minimize 2-norm(| b - A*x |)
    !using the singular value decomposition (SVD) of A. 
    !A is an M-by-N matrix which may be rank-deficient.
    svd_rcond_local = svd_rcond
    call serial_lsystem_by_svd ( phi, ymat_qr_svd, Mline, Ncolm, w_params, svd_rcond_local)

    if (allocated(ymat_qr_svd)) deallocate (ymat_qr_svd)
    if (allocated(phi)) deallocate(phi)

  case (fit_lapack_qr_constraints)

     call serial_lsystem_by_constraints(Cmat, Bmat,ymat, zmat, w_params)

  !mld_fit_type=5 ALS-Ridge block preconditioning
  case (fit_als)
    ! ALS-Ridge: works directly on the design matrix Amat (D x M)
    ! bypasses the phi = A W A^T normal equation construction.
    ! Operates on raw Amat; effective parameters are w_tilde = alpha * w
    ! (sec. 1.10.2 of MiladyNoteTechnique5.pdf).
    ! Apply sqrt(weights_snap) to Amat columns and ymat in-place so the
    ! ALS solver minimises the weighted residual  sum_j w_j ||A_j w - y_j||^2.
    do i = 1, size(Amat, 2)
      Amat(:, i) = Amat(:, i) * dsqrt(weights_snap(i))
      ymat(i, 1) = ymat(i, 1) * dsqrt(weights_snap(i))
    end do
    ! Resolve the block partition from descriptor structure (or keep user-supplied)
    call resolve_als_partition_from_descriptors(size(Amat, 1))
    block
      real(kind=kind(1.d0)), allocatable :: alpha_als(:)
      integer :: als_info, inu
      if (rangml == 0) then
        do inu = 1, als_nu_max
          call log_info("ALS-Ridge serial: block " // vtoa(inu) // " rows [" &
                        // vtoa(als_block_partition(2*inu-1)) // ", " &
                        // vtoa(als_block_partition(2*inu)) // "]")
        end do
      end if
      allocate(alpha_als(als_nu_max))
      call serial_als_ridge_solve(Amat, size(Amat,1), size(Amat,2), &
                                  ymat(:,1), als_nu_max, als_block_partition, &
                                  merge(lambda_krr, als_ridge_k, als_precond_type==als_precond_flat), &
                                  als_rho, als_nsteps, als_tol, &
                                  als_alpha_method, als_precond_type, &
                                  als_nnls_alpha, als_nnls_mode, &
                                  w_params(:,1), alpha_als, als_info)
      if (als_info /= 0 .and. rangml == 0) then
        call log_warning("ALS-Ridge serial: solver returned info = " // vtoa(als_info))
      end if
      ! Store the learned alpha_nu in the module variable for later use
      if (allocated(als_alpha)) deallocate(als_alpha)
      allocate(als_alpha(als_nu_max))
      als_alpha(:) = alpha_als(:)
      if (rangml == 0) then
        do inu = 1, als_nu_max
          call log_info("ALS-Ridge serial: alpha(" // vtoa(inu) // ") = " // vtoa(als_alpha(inu)))
        end do
      end if
      deallocate(alpha_als)
    end block

    if (allocated(ymat_qr_svd)) deallocate(ymat_qr_svd)
    if (allocated(phi)) deallocate(phi)

  case(fit_online_svd)
    ! Should not reach here — online fit has its own entry point in mld.F90
    call log_critical("fit_online_svd (mld_fit_type=6) must not use the standard serial parameter path")
    stop 'FATAL: fit_online_svd reached train_snap_get_parameters_linear'

  end select


  _MLD_END_

end subroutine train_snap_get_parameters_linear




subroutine train_snap_get_parameters_kernel()

  use ml_in_ndm_module, only: rangml,  svd_rcond, mld_fit_type, &
                              fit_home_hb, fit_lapack_qr, fit_lapack_qr_constraints, &
                              fit_lapack_ortho, fit_lapack_svd, fit_online_svd, lambda_krr
  use temporary_data_cov, only: dim_data_train
  use snap, only: w_params, Amat, ymat, Bmat, zmat, i_fit_snap, weights_snap
  use mld_logger
  use math, only: dreal_matmul, dreal_inverse
  use module_serial_linear_solver, only: serial_lsystem_by_svd, serial_lsystem_by_ortho, &
                                         serial_lsystem_by_qr, serial_lsystem_by_constraints, &
                                         serial_lsystem_by_home
  implicit none

  integer  :: i
  real(kind=kind(1.d0)), allocatable, dimension(:, :)      :: phi, phi_diag, Cmat, Cmat_transpose

  real(kind=kind(1.d0)), dimension(:, :), allocatable      :: ymat_qr_svd, Amat_big, ymat_w
  integer  :: Mline, Ncolm
  integer  :: M_obsv, P_desc
  real(kind=kind(1.d0)) :: svd_rcond_local

  _NAMECURRENT_("train_snap_get_parameters_kernel")


  _MLD_BEGIN_
  !  Amat is the transpose of the design matrix, i.e. Amat =  P_descriptor(D) x  M_observations
  !  size(Amat,1) = dim_xdesc+1, dim_xdesc**2 + dim_xdesc +1, dim_kernel etc
  !  size(Amat,2) = number of observations.
  !                      M_obsv
  !           -------------------------------
  !         |                                |
  !   Amat  |                                |   is short-fat
  !         |                                |
  !           -------------------------------

  M_obsv = size(Amat, 2)
  P_desc = size(Amat, 1)

  if (allocated(phi_diag)) deallocate (phi_diag); allocate (phi_diag(P_desc, P_desc))

  if (allocated(Cmat)) deallocate (Cmat); allocate (Cmat(P_desc, M_obsv))
  if (allocated(Cmat_transpose)) deallocate (Cmat_transpose); allocate (Cmat_transpose(M_obsv, P_desc))

  if (size(ymat, 1) /= M_obsv) then
    if (rangml == 0) write (6, *) 'Big problems in setting Amat^T w = y prroblem. The dimension of y is not the same as the number of observation on A'
    stop ' stop in train_snap_get_parameters_kernel - problem of Amat and y dimensions'
  end if

  if (i_fit_snap /= dim_data_train) then
    if (rangml == 0) write (6, *) 'The problem is not correct sized', i_fit_snap, dim_data_train
    stop 'incorrect sized problem train_snap_get_parameters'
  end if

  ! This is huge
  Cmat(:, :) = Amat(:, :)
  ! the solution is (A^T w A)^-1 A^T W ymat or (Amat w Amat^T)^-1 Amat w ymat
  ! Cmat = Amat x w
  if (allocated(ymat_w)) deallocate (ymat_w); allocate (ymat_w(size(ymat, 1), size(ymat, 2)))

  do i = 1, size(Cmat, 2)
    Cmat(:, i) = Amat(:, i)*dsqrt(weights_snap(i))
    ymat_w(i, :) = ymat(i, :)*dsqrt(weights_snap(i))
  end do
  ! This is huge


  if (lambda_krr < 0) lambda_krr = 0.d0

  ! this_take_long_time_replace_with_dgemm: phi = matmul (Amat, transpose(Cmat))
  Cmat_transpose(:, :) = transpose(Cmat)

  ! solving with regularization
  phi_diag(:, :) = 0.d0
  do i = 1, P_desc
    phi_diag(i, i) = 1.d0
  end do


  ! Premare Amat_big
  if (lambda_krr < 0) then
    ! Amat_big is just the transpose of Amat with weigths. It is the real design matrix.
    if (allocated(Amat_big)) deallocate (Amat_big); allocate (Amat_big(size(Amat, 2), size(Amat, 1)))
    Amat_big(:, :) = Cmat_transpose(:,:)
  else

    !                           Amat_big     ymat_qr_svd
    !                       (   Amat^T  )     ( y )
    !Finding w that satisfy (           ) w - (   ) =  0 is equivalent to minimize | Amat w -y |^2 + lambda |w|^2
    !                       ( lambda I  )     ( 0 )
    if (allocated(Amat_big)) deallocate (Amat_big); allocate (Amat_big(M_obsv + P_desc, P_desc))
    Amat_big(1:M_obsv, 1:P_desc) = Cmat_transpose(1:M_obsv, 1:P_desc)
    Amat_big(M_obsv + 1:M_obsv + P_desc, 1:P_desc) = lambda_krr*phi_diag(1:P_desc, 1:P_desc)
  end if


  ! Prepare ymat_qr_svd
  !  size(ymat,1) -> M_obsv
  !  size(ymat,2) -> RHS
  !  size(Amat_big,1) -> M_obsv + P_desc
  !  size(Amat_big,2) -> P_desc

  if (lambda_krr < 0) then
    if (allocated(ymat_qr_svd)) deallocate (ymat_qr_svd); allocate (ymat_qr_svd(size(ymat, 1), size(ymat, 2)))
    !k ymat_qr_svd(:, :) = ymat(:, :)
    ymat_qr_svd(:, :) = ymat_w(:, :)
  else
    if (allocated(ymat_qr_svd)) deallocate (ymat_qr_svd); allocate (ymat_qr_svd(size(ymat, 1) + P_desc, size(ymat, 2)))
    !ymat_qr_svd(1:M_obsv, :) = ymat(1:M_obsv, :)
    ymat_qr_svd(1:M_obsv, :) = ymat_w(1:M_obsv, :)
    ymat_qr_svd(1 + M_obsv:M_obsv + P_desc, :) = 0.d0
  end if

  ! Prepare lapack / scalapack
  Mline = size(Amat_big, 1)                        ! M_obsv + P_desc
  Ncolm = size(Amat_big, 2)                        ! P_desc
  !LDA = max(Mline, 1)
  !LDB = max(Mline, max(Ncolm, 1))


  select case (mld_fit_type)

  case (fit_home_hb)

    if (allocated(phi)) deallocate (phi); allocate (phi(size(Amat, 1), size(Amat, 1)))
    call dreal_matmul(Amat, size(Amat, 1), size(Amat, 2), Cmat_transpose, size(Cmat_transpose, 1), size(Cmat_transpose, 2), &
                      phi, size(phi, 1), size(phi, 2))
    phi = phi + lambda_krr*phi_diag
    ! this solve phi x w_params = Cmat x ymat -> w_params  
    !call serial_lsystem_by_home(phi, Cmat, ymat, w_params)
    call serial_lsystem_by_home(phi, ymat, Mline, Ncolm, w_params)
    deallocate (phi, Cmat, Cmat_transpose)

  case (fit_lapack_qr)

    ! QR based on Lapack dgels LLS general A MxN matrix.
    ! It is assumed that A has full rank.
    call serial_lsystem_by_qr(Amat_big, ymat_qr_svd, Mline, Ncolm, w_params)

    if (allocated(ymat_qr_svd)) deallocate (ymat_qr_svd)
    if (allocated(Amat_big)) deallocate (Amat_big)


  case (fit_lapack_ortho) 
    ! mld_fit_type=3 ortho decomposition, rank estimation
    ! ortho solution based on dgelsy
    ! DGELSY computes the minimum-norm solution to a real linear least
    ! squares problem:
    ! minimize || A * X - B ||
    ! using a complete orthogonal factorization of A.  A is an M-by-N
    ! matrix which may be rank-deficient.

    call serial_lsystem_by_ortho(Amat_big, ymat_qr_svd, Mline, Ncolm, w_params)

    if (allocated (ymat_qr_svd)) deallocate(ymat_qr_svd)
    if (allocated(Amat_big)) deallocate(Amat_big)

  case (fit_lapack_svd)   
    ! mld_fit_type=4 SVD, rank estimation
    ! svd solution based on dgelsd
    !DGELSD computes the minimum-norm solution to a real linear least
    !squares problem: minimize 2-norm(| b - A*x |)
    !using the singular value decomposition (SVD) of A. 
    !A is an M-by-N matrix which may be rank-deficient.
    svd_rcond_local = svd_rcond
    call serial_lsystem_by_svd(Amat_big, ymat_qr_svd, Mline, Ncolm, w_params, svd_rcond_local)

    if(allocated(ymat_qr_svd))  deallocate(ymat_qr_svd)
    if (allocated(Amat_big)) deallocate (Amat_big)

  case (fit_lapack_qr_constraints)

     call serial_lsystem_by_constraints(Cmat, Bmat,ymat, zmat, w_params)

  case(fit_online_svd)
    ! Should not reach here — online fit has its own entry point in mld.F90
    call log_critical("fit_online_svd (mld_fit_type=6) must not use the standard kernel parameter path")
    stop 'FATAL: fit_online_svd reached train_snap_get_parameters_kernel'

  end select

  _MLD_END_

end subroutine train_snap_get_parameters_kernel



!--------subroutines related to optimization using loss function ... no diag.

subroutine init_objective_minimization()

  use temporary_data_cov, only: dim_xdesc
  use module_objective_nl, only: lambda_nl, &
                                 jobj_T_d, jobj_E_d, jobj_F_d, jobj_S_d, &
                                 low_alpha_nl, high_alpha_nl
  use module_nlinear, only: order_nlinear, alpha_nl
  implicit none


  if (allocated(alpha_nl)) deallocate (alpha_nl); allocate (alpha_nl(order_nlinear, dim_xdesc))
  if (allocated(low_alpha_nl)) deallocate (low_alpha_nl); allocate (low_alpha_nl(order_nlinear, dim_xdesc))
  if (allocated(high_alpha_nl)) deallocate (high_alpha_nl); allocate (high_alpha_nl(order_nlinear, dim_xdesc))
  if (allocated(lambda_nl)) deallocate (lambda_nl); allocate (lambda_nl(order_nlinear))

  if (allocated(jobj_E_d)) deallocate (jobj_E_d); allocate (jobj_E_d(order_nlinear, dim_xdesc))
  if (allocated(jobj_F_d)) deallocate (jobj_F_d); allocate (jobj_F_d(order_nlinear, dim_xdesc))
  if (allocated(jobj_S_d)) deallocate (jobj_S_d); allocate (jobj_S_d(order_nlinear, dim_xdesc))

  if (allocated(jobj_T_d)) deallocate (jobj_T_d); allocate (jobj_T_d(order_nlinear*dim_xdesc))

end subroutine init_objective_minimization




subroutine train_snap_get_parameters_quadratic_objective_minimization()

  use snap, only: w_params, Amat, ymat, fit_snap
  use temporary_data_cov, only: dim_xdesc
  use derived_types, only: config_real
  use module_objective_nl, only: lambda_nl, low_alpha_nl, high_alpha_nl
  use module_nlinear, only: lnlinear_precond, order_nlinear, alpha_nl
  use module_evaluate_parameters, only: product_w_params_Amat
  implicit none
  real(kind(0.d0)), dimension(:, :), allocatable     :: Lmat
  real(kind(0.d0)), dimension(:, :), allocatable     :: ymat_copy
  real(kind(0.d0)), dimension(:, :), allocatable     :: w_params_final
  real(kind(0.d0))     :: term
  integer  :: it, nounit1, nounit2, nounit3

  if (allocated(w_params_final)) deallocate (w_params_final); allocate (w_params_final(1 + dim_xdesc + order_nlinear*dim_xdesc, 1))
  if (allocated(w_params)) deallocate (w_params); allocate (w_params(size(Amat, 1), 1))

  if (lnlinear_precond) then
    if (allocated(Lmat)) deallocate (Lmat); allocate (Lmat(dim_xdesc + 1, size(Amat, 2)))
    Lmat(1:dim_xdesc + 1, :) = Amat(1:dim_xdesc + 1, :)


    !dim_xdesc_quadratic = dim_xdesc + 1
    call train_snap_get_parameters_linear
    if (allocated(ymat_copy)) deallocate (ymat_copy); allocate (ymat_copy(size(Amat, 2), 1))
    w_params_final(1:dim_xdesc + 1, 1) = w_params(1:dim_xdesc + 1, 1)

    !Y should be copied ...
    ymat_copy = ymat

    open (file="energy_deviation_linear.milady", newunit=nounit1, action='write')
    open (file="force_deviation_linear.milady", newunit=nounit2, action='write')
    open (file="stress_deviation_linear.milady", newunit=nounit3, action='write')
    do it = 1, size(Amat, 2)
      !term = DOT_PRODUCT(w_params(:, 1), Amat(:, it))
      call product_w_params_Amat(term, it)
      if (fit_snap(it)%energy) write (nounit1, '(i9," ",3e20.10, " ", (a), " ", (a))') it, ymat(it, 1) - term, term, ymat(it, 1), config_real(fit_snap(it)%iconf)%class, trim(config_real(fit_snap(it)%iconf)%filename)
      if (fit_snap(it)%force) write (nounit2, '(i9," ",3e20.10, " ", (a), " ", (a))') it, ymat(it, 1) - term, term, ymat(it, 1), config_real(fit_snap(it)%iconf)%class, trim(config_real(fit_snap(it)%iconf)%filename)
      if (fit_snap(it)%stress) write (nounit3, '(i9," ",3e20.10, " ", (a), " ", (a))') it, ymat(it, 1) - term, term, ymat(it, 1), config_real(fit_snap(it)%iconf)%class, trim(config_real(fit_snap(it)%iconf)%filename)
      ymat(it, 1) = ymat(it, 1) - term
    end do

    close (nounit1, status='keep')
    close (nounit2, status='keep')
    close (nounit3, status='keep')


    deallocate (Amat); allocate (Amat(size(Lmat, 1), size(Lmat, 2)))
    Amat(:, :) = Lmat(:, :)

    if (allocated(w_params)) deallocate (w_params); allocate (w_params(size(Amat, 1), 1))
    call train_snap_get_parameters_linear
    low_alpha_nl(:, :) = minval(w_params)/10.d0
    high_alpha_nl(:, :) = maxval(w_params)*10.0
    ! init LBFGS ...
    alpha_nl(:, :) = 1.d-1
    lambda_nl(:) = 1.d-5
    ! start LBFGS
    call driver_lbfgs(order_nlinear*dim_xdesc)
    write (*, *) alpha_nl(1, :)
    write (*, *) alpha_nl(2, :)

    stop
  else

    w_params_final(:, :) = 0.d0
    low_alpha_nl(:, :) = -1000.d0
    high_alpha_nl(:, :) = 1000.d0

  end if

end subroutine train_snap_get_parameters_quadratic_objective_minimization



subroutine evaluate_objective()
  ! returns:
  !   -  objective function: jobj_T
  !   -  and its drivatives: jobj_T_d(:)
  !  On entry is expected that the Amat (DxM) matrix is allready filled.
  !  D is the descriptor dimension dimension:  dim_descriptor + 1
  !  M is the number of atomic environements
  !--------------------------------------------------------------
  ! There is an implementation for 2-linear
  ! In work and pending:   3-linear,  4-linear  and 5-linear
  !--------------------------------------------------------------

  use snap, only: Amat, fit_snap, ymat
  use temporary_data_cov, only: dim_xdesc
  use module_objective_nl, only: lambda_nl, n_E, n_F, n_S, &
                                 jobj_T, jobj_E, jobj_F, jobj_S, &
                                 jobj_T_d, jobj_E_d, jobj_F_d, jobj_S_d

  use module_nlinear, only: order_nlinear, alpha_nl
  
  implicit none
  real(kind(0.d0)), dimension(:, :), allocatable     :: dd_tmp_d
  real(kind(0.d0)), dimension(:), allocatable  :: dd_tmp_v
  integer  :: ii, id1, id2
  real(kind(0.d0))     :: dd_tmp
  


  n_E = 0
  n_F = 0
  n_S = 0

  jobj_T = 0.d0
  jobj_E = 0.d0
  jobj_F = 0.d0
  jobj_S = 0.d0


  jobj_T_d(:) = 0.d0
  jobj_E_d(:, :) = 0.d0
  jobj_F_d(:, :) = 0.d0
  jobj_S_d(:, :) = 0.d0

  if (order_nlinear == 2) then
    if (allocated(dd_tmp_d)) deallocate (dd_tmp_d); allocate (dd_tmp_d(order_nlinear, dim_xdesc))
    if (allocated(dd_tmp_v)) deallocate (dd_tmp_v); allocate (dd_tmp_v(order_nlinear))

    dd_tmp_d(1, 1:dim_xdesc) = 0
    dd_tmp_d(2, 1:dim_xdesc) = 0

    do ii = 1, size(Amat, 2)                         ! DxM in milady
      ! the function J for m^th column in Amat matrix.
      dd_tmp = 0.d0
      do id1 = 1, dim_xdesc
        do id2 = 1, dim_xdesc
          dd_tmp = dd_tmp + alpha_nl(1, id1)*alpha_nl(2, id2)*Amat(1 + id1, ii)*Amat(1 + id2, ii)
        end do
      end do

      ! the derivatives of J
      dd_tmp_v(1) = 2.d0*DOT_PRODUCT(alpha_nl(2, 1:dim_xdesc), Amat(2:dim_xdesc + 1, ii))
      dd_tmp_v(2) = 2.d0*DOT_PRODUCT(alpha_nl(1, 1:dim_xdesc), Amat(2:dim_xdesc + 1, ii))


      do id1 = 1, dim_xdesc
        dd_tmp_d(1, id1) = Amat(1 + id1, ii)*dd_tmp_v(1)
        dd_tmp_d(2, id1) = Amat(1 + id1, ii)*dd_tmp_v(2)
      end do


      if (fit_snap(ii)%energy) then
        n_E = n_E + 1
        jobj_E = jobj_E + fit_snap(ii)%weight*(dd_tmp - ymat(ii, 1))**2
        jobj_E_d(1, :) = jobj_E_d(1, :) + fit_snap(ii)%weight*(dd_tmp - ymat(ii, 1))*dd_tmp_d(1, :)
        jobj_E_d(2, :) = jobj_E_d(2, :) + fit_snap(ii)%weight*(dd_tmp - ymat(ii, 1))*dd_tmp_d(2, :)
      end if

      if (fit_snap(ii)%force) then
        n_F = n_F + 1
        jobj_F = jobj_F + fit_snap(ii)%weight*(dd_tmp - ymat(ii, 1))**2
        jobj_F_d(1, :) = jobj_F_d(1, :) + fit_snap(ii)%weight*(dd_tmp - ymat(ii, 1))*dd_tmp_d(1, :)
        jobj_F_d(2, :) = jobj_F_d(2, :) + fit_snap(ii)%weight*(dd_tmp - ymat(ii, 1))*dd_tmp_d(2, :)
      end if

      if (fit_snap(ii)%stress) then
        n_S = n_S + 1
        jobj_S = jobj_S + fit_snap(ii)%weight*(dd_tmp - ymat(ii, 1))**2
        jobj_S_d(1, :) = jobj_S_d(1, :) + fit_snap(ii)%weight*(dd_tmp - ymat(ii, 1))*dd_tmp_d(1, :)
        jobj_S_d(2, :) = jobj_S_d(2, :) + fit_snap(ii)%weight*(dd_tmp - ymat(ii, 1))*dd_tmp_d(2, :)
      end if


    end do                  ! im - rolled all the design matrix.

    jobj_T = 0.d0
    if (n_E > 0) then
      jobj_T = jobj_T + jobj_E/dble(n_E)
      do id1 = 1, dim_xdesc
        jobj_T_d(id1) = jobj_T_d(id1) + jobj_E_d(1, id1)/dble(n_E)
        jobj_T_d(id1 + dim_xdesc) = jobj_T_d(id1 + dim_xdesc) + jobj_E_d(2, id1)/dble(n_E)
      end do
    end if

    if (n_F > 0) then
      jobj_T = jobj_T + jobj_F/dble(n_F)
      do id1 = 1, dim_xdesc
        jobj_T_d(id1) = jobj_T_d(id1) + jobj_F_d(1, id1)/dble(n_F)
        jobj_T_d(id1 + dim_xdesc) = jobj_T_d(id1 + dim_xdesc) + jobj_F_d(2, id1)/dble(n_F)
      end do
    end if

    if (n_S > 0) then
      jobj_T = jobj_T + jobj_S/dble(n_S)
      do id1 = 1, dim_xdesc
        jobj_T_d(id1) = jobj_T_d(id1) + jobj_S_d(1, id1)/dble(n_S)
        jobj_T_d(id1 + dim_xdesc) = jobj_T_d(id1 + dim_xdesc) + jobj_S_d(2, id1)/dble(n_S)
      end do
    end if


    jobj_T = jobj_T + lambda_nl(1)*SUM(alpha_nl(1, :)**2) + lambda_nl(2)*SUM(alpha_nl(2, :)**2)
    do id1 = 1, dim_xdesc
      jobj_T_d(id1) = jobj_T_d(id1) + 2.d0*lambda_nl(1)*alpha_nl(1, id1)
      jobj_T_d(id1 + dim_xdesc) = jobj_T_d(id1 + dim_xdesc) + 2.d0*lambda_nl(2)*alpha_nl(2, id1)
    end do

  end if                  ! order_nlinear==2

end subroutine evaluate_objective




subroutine train_snap_get_parameters_quadratic_bilinear_inversion_matrix()

  use ml_in_ndm_module, only: rangml
  use temporary_data_cov, only: dim_xdesc, dim_xdesc_patch
  use derived_types, only: config_real
  use snap, only: Amat, w_params, ymat, fit_snap
  use module_mld_quadratic, only: Lmat, Qmat, dim_xdesc_quadratic, ymat_copy, w_params_final
  use module_evaluate_parameters, only: product_w_params_Amat
  use mld_logger

  implicit none

  integer  :: id1, id2, idd, lwork_svd, info_svd
  integer  :: it, size_Aline, nounit92, nounit93, nounit94
  real(kind(1.d0)), dimension(:, :), allocatable     :: w_matrix, u_svd, vt_svd
  real(kind(1.d0)), dimension(:), allocatable  :: sigma_svd, work_svd, alpha_1, alpha_2
  real(kind(1.d0))     :: term, delta_corr
  integer, parameter   :: lwmax = 1000

  _NAMECURRENT_("train_snap_get_parameters_quadratic_bilinear_inversion_matrix")


  _MLD_BEGIN_
  ! The version E - E_LML
  size_Aline = size(Amat, 1)
  call split_Amat_in_Lmat_Qmat(dim_xdesc, dim_xdesc_patch, dim_xdesc_quadratic)

  !Y should be copied ...
  if (allocated(ymat_copy)) deallocate (ymat_copy); allocate (ymat_copy(size(Amat, 2), 1))
  ymat_copy = ymat

  ! this subroutine get the parameters and ymat with the diff E_DFT - E_LML
  call precond_linear_from_Lmat

  if (allocated(w_params_final)) deallocate (w_params_final); allocate (w_params_final(size_Aline, 1))
  w_params_final(1:dim_xdesc + 1, 1) = w_params(1:dim_xdesc + 1, 1)

  deallocate (Amat); allocate (Amat(size(Qmat, 1), size(Qmat, 2)))
  Amat(:, :) = Qmat(:, :)
  if (allocated(w_params)) deallocate (w_params); allocate (w_params(size(Amat, 1), 1))

  if (allocated(w_matrix)) deallocate (w_matrix); allocate (w_matrix(dim_xdesc, dim_xdesc))
  do id1 = 1, dim_xdesc
    do id2 = 1, dim_xdesc
      idd = dim_xdesc*(id1 - 1) + id2
      w_matrix(id1, id2) = w_params(1 + dim_xdesc + idd, 1)
    end do
  end do
  !
  if (allocated(sigma_svd)) deallocate (sigma_svd); allocate (sigma_svd(size(w_matrix, 2)))
  if (allocated(u_svd)) deallocate (u_svd); allocate (u_svd(size(w_matrix, 1), size(w_matrix, 1)))
  if (allocated(vt_svd)) deallocate (vt_svd); allocate (vt_svd(size(w_matrix, 2), size(w_matrix, 2)))
  !integer   :: M-line-rows ,N - column

  !call SGESVD( JOBU, JOBVT,   M, N, A,  LDA, S,U, LDU, VT, LDVT, WORK, LWORK, INFO )
  if (allocated(work_svd)) deallocate (work_svd); allocate (work_svd(lwmax))
  lwork_svd = -1
  call dgesvd('A', 'A', size(w_matrix, 1), size(w_matrix, 2), w_matrix, size(w_matrix, 1), &
              sigma_svd, u_svd, size(w_matrix, 1), vt_svd, size(w_matrix, 2), work_svd, lwork_svd, info_svd)
  if (info_svd < 0) then
    if (rangml == 0) then
      write (6, *) '1st call - Illegal setup for SVD in parameters.F90 the info has the value ', info_svd
    end if
  end if

  lwork_svd = int(work_svd(1))
  if (allocated(work_svd)) deallocate (work_svd); allocate (work_svd(lwork_svd))
  call dgesvd('A', 'A', size(w_matrix, 1), size(w_matrix, 2), w_matrix, size(w_matrix, 1), &
              sigma_svd, u_svd, size(w_matrix, 1), vt_svd, size(w_matrix, 2), work_svd, lwork_svd, info_svd)

  if (info_svd < 0) then
    if (rangml == 0) then
      write (6, *) '2nd call - Illegal setup for SVD in parameters.F90 the info has the value ', info_svd
    end if
  end if

  if (allocated(alpha_1)) deallocate (alpha_1); allocate (alpha_1(size(w_matrix, 1)))
  if (allocated(alpha_2)) deallocate (alpha_2); allocate (alpha_2(size(w_matrix, 1)))
  alpha_1(:) = u_svd(:, 1)
  alpha_2(:) = vt_svd(1, :)*sigma_svd(1)
  do id1 = 1, dim_xdesc
    do id2 = 1, dim_xdesc
      idd = dim_xdesc*(id1 - 1) + id2
      w_params(idd, 1) = alpha_1(id1)*alpha_2(id2)
    end do
  end do
  w_params_final(dim_xdesc + 2:dim_xdesc_quadratic, 1) = w_params(1:size(Amat, 1), 1)

  open (file="energy_deviation_quadratic.milady", newunit=nounit92, action='write')
  open (file="force_deviation_quadratic.milady", newunit=nounit93, action='write')
  open (file="stress_deviation_quadratic.milady", newunit=nounit94, action='write')
  do it = 1, size(Amat, 2)
    !term = DOT_PRODUCT(w_params(:, 1), Amat(:, it))
    call product_w_params_Amat(term, it)
    delta_corr = ymat(it, 1) - term
    if (fit_snap(it)%energy) write (nounit92, '(i9," ",3e20.10, " ", (a), " ", (a))') it, delta_corr, term, ymat(it, 1), config_real(fit_snap(it)%iconf)%class, trim(config_real(fit_snap(it)%iconf)%filename)
    if (fit_snap(it)%force) write (nounit93, '(i9," ",3e20.10, " ", (a), " ", (a))') it, delta_corr, term, ymat(it, 1), config_real(fit_snap(it)%iconf)%class, trim(config_real(fit_snap(it)%iconf)%filename)
    if (fit_snap(it)%stress) write (nounit94, '(i9," ",3e20.10, " ", (a), " ", (a))') it, delta_corr, term, ymat(it, 1), config_real(fit_snap(it)%iconf)%class, trim(config_real(fit_snap(it)%iconf)%filename)
  end do
  close (nounit92, status='keep')
  close (nounit93, status='keep')
  close (nounit94, status='keep')

  ymat = ymat_copy
  deallocate (Amat); allocate (Amat(dim_xdesc_quadratic, size(Qmat, 2)))
  Amat(1:dim_xdesc + 1, :) = Lmat(1:dim_xdesc + 1, :)
  deallocate (Lmat)
  Amat(dim_xdesc + 2:dim_xdesc_quadratic, :) = Qmat(:, :)
  deallocate (Qmat)
  if (allocated(w_params)) deallocate (w_params); allocate (w_params(dim_xdesc_quadratic, 1))
  w_params(1:dim_xdesc_quadratic, :) = w_params_final(1:dim_xdesc_quadratic, :)
  deallocate (w_params_final)

  _MLD_END_

end subroutine train_snap_get_parameters_quadratic_bilinear_inversion_matrix


!-------------------------------------------------------------------------------
! resolve_als_partition_from_descriptors
!
! Builds the block partition for ALS-Ridge from the descriptor structure.
! If the user already supplied a non-zero partition (via some external
! mechanism), it is left untouched.  Otherwise the partition is inferred from
! the current descriptor_type:
!
!   descriptor_body  —  one block per active body order (l_body_order == .true.)
!                       block dimension = dim_desc_body(i)
!
!   descriptor_ace   —  one block per active ACE body order (l_ace_order == .true.)
!                       body 1: base_params(1)%mumax * (nkmax_order1+1-nkmin_order1)
!                       body n>=2: size(base_cnlm(n)%dicoB) * base_params(n)%kmax
!                       each block is replicated zetaace_order times when > 1.
!
!   anything else    —  single block [1, D_total]
!
! Arguments:
!   D_total  — total number of descriptor rows (== dim_xdesc, dimr_sca_Amat, or
!              size(Amat,1) depending on the calling path).
!-------------------------------------------------------------------------------
subroutine resolve_als_partition_from_descriptors(D_total)

#include "../MLD_MACROS.INC"

  use ml_in_ndm_module, only: descriptor_type, descriptor_body, descriptor_ace, descriptor_ftnbody
  use module_input_als_fit, only: als_block_partition, als_nu_max, &
                                  set_als_block_partition
  use mld_logger

  implicit none
  integer, intent(in) :: D_total

  _NAMECURRENT_("resolve_als_partition_from_descriptors")

  _MLD_BEGIN_

  select case (descriptor_type)

    case (descriptor_body)
      call resolve_als_partition_body(D_total)

    case (descriptor_ace)
      call resolve_als_partition_ace(D_total)

    case (descriptor_ftnbody)
      call resolve_als_partition_ftnbody(D_total)

    case default
      ! Other descriptor types: single-block default
      call set_als_block_partition(D_total, 0)

  end select

  _MLD_END_

end subroutine resolve_als_partition_from_descriptors


!-------------------------------------------------------------------------------
! resolve_als_partition_body
!
! Build ALS block partition for k2b (body) descriptors.
! One block per active body order using dim_desc_body(1:6) / l_body_order(1:6).
!-------------------------------------------------------------------------------
subroutine resolve_als_partition_body(D_total)

#include "../MLD_MACROS.INC"

  use module_body_desc, only: dim_desc_body, l_body_order
  use module_input_als_fit, only: set_als_block_partition
  use mld_logger

  implicit none
  integer, intent(in) :: D_total

  integer :: nblocks, i
  integer, allocatable :: block_dims(:)

  _NAMECURRENT_("resolve_als_partition_body")

  _MLD_BEGIN_

  ! Count active body orders with non-zero dimension
  nblocks = 0
  do i = 1, 6
    if (l_body_order(i) .and. dim_desc_body(i) > 0) nblocks = nblocks + 1
  end do

  if (nblocks == 0) then
    ! No active body orders — fall back to single block
    call set_als_block_partition(D_total, 0)
    _MLD_END_
    return
  end if

  ! +1 for the bias row (row 1 of Amat = nat)
  allocate(block_dims(nblocks + 1))
  block_dims(1) = 1
  nblocks = 0
  do i = 1, 6
    if (l_body_order(i) .and. dim_desc_body(i) > 0) then
      nblocks = nblocks + 1
      block_dims(nblocks + 1) = dim_desc_body(i)
    end if
  end do
  nblocks = nblocks + 1

  call log_info("ALS-Ridge: body descriptor partition with " // vtoa(nblocks) // " blocks")

  call set_als_block_partition(D_total, nblocks, block_dims)
  deallocate(block_dims)

  _MLD_END_

end subroutine resolve_als_partition_body


!-------------------------------------------------------------------------------
! resolve_als_partition_ace
!
! Build ALS block partition for ACE descriptors.
! One block per active ACE body order (times zetaace_order replicas).
!
! Memory layout of the ACE descriptor in the design matrix rows:
!
!   ACE_CHEM_INCOMPLETE (single species):
!     [ body1_z1 | body1_z2 | ... | body2_z1 | body2_z2 | ... ]
!     Total = ace_dim = SUM over active body orders of (dim_i * n_zeta)
!
!   ACE_CHEM_STANDARD (multi-species):
!     [ species1: body1_z1|..|bodyN_zZ | species2: body1_z1|..|bodyN_zZ | ... ]
!     Total = ace_dim * fix_no_of_elements
!     Each species block has the same internal body-order/zeta structure.
!
! For the multi-species case the partition lists all body-order×zeta blocks
! of species 1 first, then species 2, etc., giving one contiguous block per
! (species, body_order, zeta) combination.
!
! Per body order:
!   body 1  : dim = mumax * (nkmax_order1 + 1 - nkmin_order1)
!   body n>=2: dim = size(base_cnlm(n)%dicoB) * kmax
!-------------------------------------------------------------------------------
subroutine resolve_als_partition_ace(D_total)

#include "../MLD_MACROS.INC"

  use module_ace_desc, only: ace_numax, l_ace_order, zetaace_order, &
                             ace_radial_chem, nkmax_order1, nkmin_order1, &
                             base_params, ace_chem, &
                             ACE_CHEM_INCOMPLETE, ACE_CHEM_STANDARD
  use module_base_cnlm, only: base_cnlm
  use module_chemical_species, only: fix_no_of_elements
  use module_input_als_fit, only: set_als_block_partition
  use mld_logger

  implicit none
  integer, intent(in) :: D_total

  integer :: nactive, i, iz, isp, bdim, idx, n_zeta, n_species, total_blocks
  integer, allocatable :: block_dims(:)

  _NAMECURRENT_("resolve_als_partition_ace")

  _MLD_BEGIN_

  n_zeta = max(1, zetaace_order)

  ! Determine species multiplier
  if (ace_chem == ACE_CHEM_STANDARD) then
    n_species = fix_no_of_elements
  else
    n_species = 1
  end if

  ! Count active body orders
  nactive = 0
  do i = 1, ace_numax
    if (l_ace_order(i)) nactive = nactive + 1
  end do

  if (nactive == 0) then
    call set_als_block_partition(D_total, 0)
    _MLD_END_
    return
  end if

  ! Total blocks = 1 (bias row) + n_species * nactive * n_zeta
  ! Row 1 of Amat is the bias term (nat), so we add a leading block of size 1.
  total_blocks = 1 + n_species * nactive * n_zeta
  allocate(block_dims(total_blocks))

  ! Block 1: bias row (row 1 of Amat = nat)
  idx = 1
  block_dims(1) = 1

  ! The layout is: [bias] [species1: body_orders x zeta] [species2: ...] ...
  do isp = 1, n_species
    do i = 1, ace_numax
      if (.not. l_ace_order(i)) cycle
      ! Compute the B-basis dimension for this body order
      if (i == 1) then
        bdim = base_params(1)%mumax * (nkmax_order1 + 1 - nkmin_order1)
      else
        bdim = size(base_cnlm(i)%dicoB) * base_params(i)%kmax
      end if
      ! Replicate for each zeta power
      do iz = 1, n_zeta
        idx = idx + 1
        block_dims(idx) = bdim
      end do
    end do
  end do

  call log_info("ALS-Ridge: ACE descriptor partition with " // vtoa(idx) &
                // " blocks (nu_max=" // vtoa(ace_numax) &
                // ", zeta_order=" // vtoa(n_zeta) &
                // ", n_species=" // vtoa(n_species) // ")")

  call set_als_block_partition(D_total, idx, block_dims)
  deallocate(block_dims)

  _MLD_END_

end subroutine resolve_als_partition_ace

!-------------------------------------------------------------------------------
! resolve_als_partition_ftnbody
!
! Build the ALS-Ridge block partition for the (multispecies) ftnbody descriptor.
! The descriptor is naturally block-structured: per active body order n, a stack
! of ftnbody_n_channels(n) chemical channels, each dim_rff(n) wide
! (dim_desc_body(n) = dim_rff(n) * ftnbody_n_channels(n)).  We expose a
! mode-aware, contiguous block partition so ALS-Ridge can renormalise / ridge
! each chemical group on its own scale:
!   bias row                                   -> 1 block
!   mode 0 (single species)                    -> 1 block per active order
!   mode 1 (exact)   group by central species  -> S blocks per order  (each = dim_desc_body/S)
!   mode 2 (low-rank)                          -> R_n blocks per order (each = dim_rff(n))
!   mode 3 (hash)                              -> H_n blocks per order (each = dim_rff(n))
! Grouping the exact channels by central species keeps the block count = S (not
! N_chem) and avoids empty blocks on a dataset that contains every species as a
! central atom.  The solver guards any residual empty block (alpha=1, lambda=0).
!-------------------------------------------------------------------------------
subroutine resolve_als_partition_ftnbody(D_total)

#include "../MLD_MACROS.INC"

  use module_body_desc, only: dim_desc_body, l_body_order
  use module_ftnbody,   only: dim_rff, ftnbody_chem_mode, ftnbody_n_channels
  use module_chemical_species, only: fix_no_of_elements
  use module_input_als_fit, only: set_als_block_partition
  use mld_logger

  implicit none
  integer, intent(in) :: D_total

  integer :: n, b, nblk_n, blksz_n, nblocks, idx
  integer, allocatable :: block_dims(:)

  _NAMECURRENT_("resolve_als_partition_ftnbody")

  _MLD_BEGIN_

  ! Count blocks: 1 (bias) + per active order its channel groups.
  nblocks = 1
  do n = 2, 5
    if (l_body_order(n) .and. dim_desc_body(n) > 0) then
      call ftnbody_block_layout(n, nblk_n, blksz_n)
      nblocks = nblocks + nblk_n
    end if
  end do

  if (nblocks <= 1) then
    ! nothing active -> single-block default
    call set_als_block_partition(D_total, 0)
    _MLD_END_
    return
  end if

  allocate(block_dims(nblocks))
  block_dims(1) = 1        ! bias row (row 1 of Amat = nat)
  idx = 1
  do n = 2, 5
    if (l_body_order(n) .and. dim_desc_body(n) > 0) then
      call ftnbody_block_layout(n, nblk_n, blksz_n)
      do b = 1, nblk_n
        idx = idx + 1
        block_dims(idx) = blksz_n
      end do
    end if
  end do

  call log_info("ALS-Ridge: ftnbody chem-block partition, mode " // vtoa(ftnbody_chem_mode) &
                // ", " // vtoa(nblocks) // " blocks")

  call set_als_block_partition(D_total, nblocks, block_dims)
  deallocate(block_dims)

  _MLD_END_

contains

  ! number (nblk) and width (blksz) of the ALS blocks for one active body order
  subroutine ftnbody_block_layout(nn, nblk, blksz)
    integer, intent(in)  :: nn
    integer, intent(out) :: nblk, blksz
    select case (ftnbody_chem_mode)
    case (0)            ! single species: one block per order
      nblk  = 1
      blksz = dim_desc_body(nn)
    case (1)            ! exact: group contiguous channels by central species
      nblk  = fix_no_of_elements
      blksz = dim_desc_body(nn) / fix_no_of_elements
    case default        ! low-rank / hash: one block per channel (rank or bucket)
      nblk  = ftnbody_n_channels(nn)
      blksz = dim_rff(nn)
    end select
  end subroutine ftnbody_block_layout

end subroutine resolve_als_partition_ftnbody
