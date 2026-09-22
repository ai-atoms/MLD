! HND XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
! HND X
! HND X   MiLaDy (Machine Learning Dynamics) 
! HND X   
! HND X
! HND X   Milady was written and designed by Mihai-Cosmin Marinica, Alexandra M. Goryaeva 
! HND X   Copyright 2015-2023.
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

module module_objective_nl

  real(kind(0.d0)), dimension(:), allocatable  :: lambda_nl
  real(kind(0.d0))     :: jobj_T, jobj_E, jobj_F, jobj_S
  real(kind(0.d0)), dimension(:), allocatable  :: jobj_T_d
  real(kind(0.d0)), dimension(:, :), allocatable     :: jobj_E_d, jobj_F_d, jobj_S_d
  integer  :: n_E, n_F, n_S
  real(kind(0.d0)), dimension(:, :), allocatable     :: low_alpha_nl, high_alpha_nl

end module module_objective_nl


module snap_interface
  !$-------------------------------------------------------------
  ! Interfaces for some functions used in snap
  ! TODO TO REMOVE
  !
  !
  !
  !
  !$-------------------------------------------------------------
  interface

    subroutine train_fill_Amat_with_force(iconf, pack_opt)
      integer, intent(in)  :: iconf
      logical, intent(in), optional    :: pack_opt
    end subroutine

    subroutine md_mld_compute_energy(iconf, pack_opt)
      integer, intent(in)  :: iconf
      logical, intent(in), optional    :: pack_opt
    end subroutine

    subroutine md_mld_compute_force(iconf, pack_opt)
      integer, intent(in)  :: iconf
      logical, intent(in), optional    :: pack_opt
    end subroutine

    subroutine md_mld_compute_stress(iconf, pack_opt)
      integer, intent(in)  :: iconf
      logical, intent(in), optional    :: pack_opt
    end subroutine

  end interface

end module snap_interface

module snap
  use module_kind_variables, only: kind_double
  ! Store the tools for snap fitting
  ! A matrix  from Amat^T x  w_params = ymat.
  ! In the case of LML:
  !       Amat    (1 + dim_xdesc , dim_data_train)
  !       ymat    (dim_data_train,              1)
  !       w_params(1 + dim_xdesc ,              1)
  ! training Amat and ymat

  integer  :: dim_xdesc_linear, dim_xdesc_full, dim_design_line, dim_design_column
  real(kind=kind(1.d0)), dimension(:, :), allocatable      :: Amat
  real(kind=kind(1.d0)), dimension(:, :), allocatable      :: ymat
  !constraints Bmat and zmat
  real(kind=kind(1.d0)), dimension(:, :), allocatable      :: Bmat
  real(kind=kind(1.d0)), dimension(:, :), allocatable      :: zmat

  real(kind=kind(1.d0)), dimension(:), allocatable   :: weights_snap
  ! Per-datapoint ZBL contribution (energy/force/stress), gathered globally
  ! like ymat. It stores the ZBL part that is subtracted from ymat during
  ! design-matrix construction, so the train report can add it back without
  ! touching the subworld-local config_real%ezbl/%fzbl/%szbl (which are not
  ! gathered). Zero unless zbl_potential .and. zbl_type == zbl_mode_alone.
  real(kind=kind(1.d0)), dimension(:), allocatable   :: y_zbl_train
  !parameters
  real(kind=kind(1.d0)), dimension(:, :), allocatable      :: w_params

  ! MD part of the code
  real(kind(0.d0))     :: ene_snap
  real(kind(0.d0)), dimension(:, :), allocatable     :: fp_snap
  real(kind(0.d0)), dimension(:, :), allocatable     :: fp_var
  real(kind(0.d0)), dimension(6)   :: stress_snap
  ! D, dimension 
  real(kind_double), dimension(:),     allocatable :: xdesc_emd 
  ! D, N, 3 dimension
  real(kind_double), dimension(:,:,:), allocatable :: xdesc_fmd
  ! D, 6 dimension
  real(kind_double), dimension(:,:), allocatable :: xdesc_smd
  ! end MD 
  integer  :: i_fit_snap, i_constraints_snap

  integer  :: dim_ene_train_lml, dim_force_train_lml, dim_stress_train_lml
  integer  :: dim_ene_constraints, dim_force_constraints, dim_stress_constraints

  integer  :: i_e_train_snap, i_f_train_snap, i_s_train_snap
  real(kind(0.d0)), dimension(:), allocatable  :: y_e_train_snap, y_f_train_snap, y_s_train_snap
  real(kind(0.d0)), dimension(:), allocatable  :: y_e_train_base, y_f_train_base, y_s_train_base
  real(kind(0.d0)), dimension(:), allocatable  :: y_e_p_a_train_base, y_e_p_a_train_snap

  !testing
  integer  :: dim_ene_test_snap, dim_force_test_snap, dim_stress_test_snap
  real(kind(0.d0)), dimension(:), allocatable  :: y_e_test_snap, y_f_test_snap, y_s_test_snap, var_f_test_snap
  real(kind(0.d0)), dimension(:), allocatable  :: y_e_test_base, y_f_test_base, y_s_test_base
  real(kind(0.d0)), dimension(:), allocatable  :: y_e_p_a_test_base, y_e_p_a_test_snap
  integer, dimension(:), allocatable :: e_test_snap, s_test_snap, f_test_snap


  real(kind(0.d0))     :: train_rmse_energy, train_rmse_force, train_rmse_stress, &
                          train_mae_energy, train_mae_force, train_mae_stress, &
                          test_rmse_energy, test_rmse_force, test_rmse_stress, &
                          test_mae_energy, test_mae_force, test_mae_stress

  ! the original weigths vector of dimension dim_data_train
  real(kind=kind(1.d0)), dimension(:), allocatable   :: tmp_weights_snap
  real(kind(0.d0))     :: Jfunc

  type fit_snap_type
    integer  :: iconf
    logical  :: force = .false.
    logical  :: energy = .false.
    logical  :: stress = .false.
    real(kind(0.d0))     :: weight
    integer  :: iatom
    integer  :: ix  
  end type fit_snap_type

  type weights_type
    integer  :: db_line
    logical  :: has_force = .false.
    logical  :: has_energy = .false.
    logical  :: has_stress = .false.
  end type weights_type

  type db_type
    integer  :: i_e
    integer  :: i_f
    integer  :: i_s
  end type db_type

  type(fit_snap_type), dimension(:), allocatable     :: fit_snap
  type(weights_type), dimension(:), allocatable      :: map_weights_in_db
  type(db_type), dimension(:), allocatable     :: map_db_in_weights

end module snap




module module_mld_quadratic
  ! module dedicated for QML/QNML implementation using the LML framework.
  ! dim_xdesc_kernel = 1 + dim_xdesc (linear_part) + dim_quadratic  ( cim_xdec**2 quadratic part)
  ! Lmat = (1 + dim_xdesc, dim_data_train = size(Amat,2)) = linear part of the design matrix
  !      --> matrix used in parameters.F90 subroutines:
  ! Qmat = (dim_quadratic, dim_data_train = size(Amat,2))    = kernel part of the design matrix  .
  ! In the case of kernel:
  !       ymat_copy:     (dim_data_train_snap = size(Amat,2), 1 ))
  !       sub_yfunc_E_train   (dim_ene_train_lml,              1)
  !       w_params_final (dim_xdesc_kernel = 1 + dim_xdesc + dim_kernel)
  !--------------------------------------------------------------
  !       Not used for the moement:
  !--------------------------------------------------------------
  !in order to make the link between the design matrix and fit.
  !       map_Emat_index_to_fit
  !       map_Fmat_index_to_fit
  !       map_Smat_index_to_fit
  !       sub_yfunc_E_train
  !       sub_yfunc_F_train
  !       sub_yfunc_S_train
  ! training Emat and ymat_E

  integer  :: dim_xdesc_quadratic, dim_quadratic
  real(kind=kind(1.d0)), dimension(:, :), allocatable      :: Lmat
  real(kind=kind(1.d0)), dimension(:, :), allocatable      :: Qmat
  real(kind=kind(1.d0)), dimension(:, :), allocatable      :: ymat_copy
  real(kind=kind(1.d0)), dimension(:, :), allocatable      :: w_params_final
  !real(kind=kind(1.d0)), dimension(:,:), allocatable :: ymat_E
  real(kind=kind(1.d0)), dimension(:), allocatable   :: sub_yfunc_E_train
  integer, dimension(:), allocatable     :: map_Emat_index_to_fit

  ! Fmat has dimension (1+dim_desc)^2 x n_data
  !real(kind=kind(1.d0)), dimension(:,:), allocatable :: ymat_F
  real(kind=kind(1.d0)), dimension(:), allocatable   :: sub_yfunc_F_train
  integer, dimension(:), allocatable     :: map_Fmat_index_to_fit


  ! Fmat has dimension D^2 x n_data
  !real(kind=kind(1.d0)), dimension(:,:), allocatable :: Smat
  !real(kind=kind(1.d0)), dimension(:,:), allocatable :: ymat_S
  real(kind=kind(1.d0)), dimension(:), allocatable   :: sub_yfunc_S_train
  integer, dimension(:), allocatable     :: map_Smat_index_to_fit


  integer  :: i_e_fit_snap, i_f_fit_snap, i_s_fit_snap
  !parameters
  !real(kind=kind(1.d0)), dimension(:,:), allocatable :: w_params

end module module_mld_quadratic



module module_mld_polyc
  ! module dedicated for PolyChaos implementation using the LML framework.
  ! dim_xdesc_polyc = 1 + dim_xdesc (linear_part) + ...
  ! Lmat = (1 + dim_xdesc, dim_data_train = size(Amat,2)) = linear part of the design matrix
  !      --> matrix used in parameters.F90 subroutines:
  ! Qmat = (dim_kernel, dim_data_train = size(Amat,2))    = kernel part of the design matrix  .
  ! In the case of kernel:
  !       ymat_copy:     (dim_data_train_snap = size(Amat,2), 1 ))
  !       sub_yfunc_E_train   (dim_ene_train_lml,              1)
  !       w_params_final (dim_xdesc_kernel = 1 + dim_xdesc + dim_kernel)
  !--------------------------------------------------------------
  !       Not used for the moment:
  !--------------------------------------------------------------
  ! in order to make the link between the design matrix and fit.
  !       map_Emat_index_to_fit
  !       map_Fmat_index_to_fit
  !       map_Smat_index_to_fit

  integer  :: dim_xdesc_polyc, dim_polyc
  real(kind=kind(1.d0)), dimension(:, :), allocatable      :: Lmat
  real(kind=kind(1.d0)), dimension(:, :), allocatable      :: Qmat
  real(kind=kind(1.d0)), dimension(:, :), allocatable      :: ymat_copy
  real(kind=kind(1.d0)), dimension(:, :), allocatable      :: w_params_final


  integer, dimension(:), allocatable     :: map_Emat_index_to_fit
  integer, dimension(:), allocatable     :: map_Fmat_index_to_fit
  integer, dimension(:), allocatable     :: map_Smat_index_to_fit


  integer  :: i_e_fit_snap, i_f_fit_snap, i_s_fit_snap

end module module_mld_polyc




module module_mld_kernel
  ! module dedicated for Kernel implementation using the LML framework.
  ! dim_xdesc_kernel = 1 + dim_xdesc (linear_part) + dim_kernel (kernel part)
  ! Lmat = (1 + dim_xdesc, dim_data_train = size(Amat,2)) = linear part of the design matrix
  !      --> matrix used in parameters.F90 subroutines:
  ! Qmat = (dim_kernel, dim_data_train = size(Amat,2))    = kernel part of the design matrix  .
  ! In the case of kernel:
  !       ymat_copy:     (dim_data_train_snap = size(Amat,2), 1 ))
  !       sub_yfunc_E_train   (dim_ene_train_lml,              1)
  !       w_params_final (dim_xdesc_kernel = 1 + dim_xdesc + dim_kernel)
  !--------------------------------------------------------------
  !       Not used for the moement:
  !--------------------------------------------------------------
  !in order to make the link between the design matrix and fit.
  !       map_Emat_index_to_fit
  !       map_Fmat_index_to_fit
  !       map_Smat_index_to_fit


  integer  :: dim_xdesc_kernel

  real(kind=kind(1.d0)), dimension(:, :), allocatable      :: Lmat
  real(kind=kind(1.d0)), dimension(:, :), allocatable      :: Qmat
  real(kind=kind(1.d0)), dimension(:, :), allocatable      :: ymat_copy
  real(kind=kind(1.d0)), dimension(:, :), allocatable      :: w_params_final
  !real(kind=kind(1.d0)), dimension(:,:), allocatable :: ymat_E
  integer, dimension(:), allocatable     :: map_Emat_index_to_fit
  integer, dimension(:), allocatable     :: map_Fmat_index_to_fit
  integer, dimension(:), allocatable     :: map_Smat_index_to_fit

  integer  :: i_e_fit_snap, i_f_fit_snap, i_s_fit_snap

end module module_mld_kernel



module module_snap_nlinear
  !       ymat_copy:     (dim_data_train_snap = size(Amat,2), 1 ))
  !       sub_yfunc_E_train   (dim_ene_train_lml,              1)
  !       w_params_final (dim_xdesc_nlinear = 1 + dim_xdesc + dim_linear)
  !--------------------------------------------------------------
  !       Not used for the moment:
  !--------------------------------------------------------------
  !in order to make the link between the design matrix and fit.
  !       map_Emat_index_to_fit
  !       map_Fmat_index_to_fit
  !       map_Smat_index_to_fit


  integer  :: dim_xdesc_nlinear, dim_nlinear

  real(kind=kind(1.d0)), dimension(:, :), allocatable      :: Lmat
  real(kind=kind(1.d0)), dimension(:, :), allocatable      :: Qmat
  real(kind=kind(1.d0)), dimension(:, :), allocatable      :: ymat_copy
  real(kind=kind(1.d0)), dimension(:, :), allocatable      :: w_params_final
  ! real(kind=kind(1.d0)), dimension(:,:), allocatable :: ymat_E
  integer, dimension(:), allocatable     :: map_Emat_index_to_fit
  integer, dimension(:), allocatable     :: map_Fmat_index_to_fit
  integer, dimension(:), allocatable     :: map_Smat_index_to_fit

  integer  :: i_e_fit_snap, i_f_fit_snap, i_s_fit_snap

end module module_snap_nlinear
