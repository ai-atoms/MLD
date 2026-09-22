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

! WEIGHTS OPTIMIZATION

#include "../MLD_MACROS.INC"



module mod_function 
  contains 
subroutine function_to_min(xval, objval)
  use data_type, only: nunit
  use ml_in_ndm_module, only: rangml, debug, lambda_krr_2, lambda_krr
  use module_optimization, only:  factor_energy_error, factor_force_error, factor_stress_error,  &
                              optimize_weights_L1, optimize_weights_L2, optimize_weights_Le
  use snap, only:  weights_snap, tmp_weights_snap
  use module_optimization, only: lambda_krr_fake, dim_weights_function, dim_weights_function_full, & 
                                 tmp_weights, itopt 
  implicit none
  real(kind(0.d0)), dimension(dim_weights_function + 1), intent(in)    :: xval
  real(kind(0.d0)), intent(out)    :: objval
  real(kind(0.d0))     :: mae_energy, mae_force, mae_stress
  tmp_weights(1:dim_weights_function) = xval(1:dim_weights_function)
  call fill_tmp_weigths_with_w

  weights_snap(:) = tmp_weights_snap(:)
  lambda_krr = lambda_krr_fake
  lambda_krr = xval(dim_weights_function_full)
  call train_snap_get_parameters
  lambda_krr = xval(dim_weights_function_full)
  if (optimize_weights_Le) lambda_krr_2 = xval(dim_weights_function_full - 1)
  call train_minimize_error(mae_energy, mae_force, mae_stress)
  !debug if (rangml==0) write (6,'("44444", 6e20.10)') factor_energy_error, mae_energy, factor_force_error, mae_force, factor_stress_error, mae_stress
  itopt = itopt + 1
  objval = factor_energy_error*mae_energy + factor_force_error*mae_force + factor_stress_error*mae_stress
  if (lambda_krr >= 0) then
    if (optimize_weights_L2) objval = objval + lambda_krr*sum(weights_snap(:)**2)
    if (optimize_weights_L1) objval = objval + lambda_krr*sum(dabs(weights_snap(:)))
    if (optimize_weights_Le) objval = objval + lambda_krr_2*sum(dabs(weights_snap(:))) + lambda_krr*sum(weights_snap(:)**2)
  end if
  if (debug) then

    if (optimize_weights_Le) then
      if (rangml == 0) write (6, '("ML:   Jfunc optimization :",i7, 3e20.10)') itopt, lambda_krr, lambda_krr_2, objval, sum(weights_snap(:)**2)
    else
      if (rangml == 0) write (6, '("ML:   Jfunc optimization :",i7, 3e20.10)') itopt, lambda_krr, objval, sum(weights_snap(:)**2)
    end if
  end if
  write (nunit, '(i7,4e20.10)') itopt, objval, mae_energy, mae_force, mae_stress
end subroutine function_to_min
end module mod_function

module hyper_optimize_chemical 
  implicit none
  
  contains

  subroutine get_dim_objective_function_chemical(dim_chemical_function)
    use ml_in_ndm_module, only: weighted, weighted_3ch
    use module_chemical_species, only: fix_no_of_elements        
    use mld_logger
    use mld_string
    use mld_mpi
  
    implicit none 
    integer, intent(inout) :: dim_chemical_function
  
    _NAMECURRENT_("get_dim_objective_function_chemical")
  
    _MLD_BEGIN_
  
    if (.not.(weighted)) then 
      call log_warning('There is nothing to optimize when the system is not weighted=.true.')
      call mld_mpi_abort('MLD FATAL: this crazy run is stoped') 
    end if 
    
    if (fix_no_of_elements <= 1) then 
      call log_warning('There is nothing to optimize for chemical weights for a single element')
      call mld_mpi_abort('MLD FATAL: this crazy run is stoped') 
    end if 
  
    dim_chemical_function = 0 

    if (weighted) then 
       dim_chemical_function = fix_no_of_elements 
    end if
  
    if (weighted_3ch) then
      dim_chemical_function = dim_chemical_function + fix_no_of_elements
    end if 
  
    return
  
    _MLD_END_
  end  subroutine get_dim_objective_function_chemical 
  
  subroutine init_from_chemical_in_objective_function (dim_chemical_function, tmp_chemical, upper_chemical, lower_chemical)
    use module_kind_variables, only: kind_double
    use ml_in_ndm_module, only: weighted, weighted_3ch, fix_weighted_for_element_ini, fix_weighted_for_element_3ch_ini
    use module_chemical_species, only: fix_no_of_elements                      
    use mld_logger

    implicit none
    integer, intent(in) :: dim_chemical_function  
    real(kind_double), dimension(:), allocatable, intent(inout) :: tmp_chemical, upper_chemical, lower_chemical
  
    _NAMECURRENT_("init_from_chemical_in_objective_function")
  
    _MLD_BEGIN_

    if (allocated(tmp_chemical)) deallocate (tmp_chemical)  
    allocate (tmp_chemical(dim_chemical_function))
    if (allocated(lower_chemical)) deallocate (lower_chemical) 
    allocate (lower_chemical(dim_chemical_function))
    if (allocated(upper_chemical)) deallocate (upper_chemical) 
    allocate (upper_chemical(dim_chemical_function))
  
    lower_chemical(:) =  0.2d0 
    upper_chemical(:) =  1.5d0

    if (weighted) then 
      tmp_chemical(1:fix_no_of_elements) = fix_weighted_for_element_ini(:)
      if (weighted_3ch) then 
        tmp_chemical(fix_no_of_elements+1:dim_chemical_function) = fix_weighted_for_element_3ch_ini
      end if   
    end if 

    _MLD_END_ 
  end subroutine init_from_chemical_in_objective_function

  subroutine function_to_min_chemical(xval, objval)
    use data_type, only: nunit
    use ml_in_ndm_module, only: rangml, weighted, weighted_3ch, &
                                fix_weighted_for_element, fix_weighted_for_element_3ch, &
                                fix_weighted_for_element_ini, fix_weighted_for_element_3ch_ini, &
                                weighted_auto, rangml 
    use module_chemical_species, only: fix_no_of_elements

    use module_optimization, only:  factor_energy_error, factor_force_error, factor_stress_error, & 
                                    dim_chemical_function, tmp_chemical, itopt 
    use module_test_mld, only: main_test_mld
    use main_mld_mod, only: main_train_mld
    use snap, only: train_mae_energy, train_mae_force, train_mae_stress                                
    implicit none
    real(kind(0.d0)), dimension(dim_chemical_function), intent(in)    :: xval
    real(kind(0.d0)), intent(out)    :: objval
    real(kind(0.d0))     :: mae_energy, mae_force, mae_stress

    tmp_chemical(1:dim_chemical_function) = xval(1:dim_chemical_function)
    weighted_auto=.false.
    if (weighted) then 
      fix_weighted_for_element_ini(1:fix_no_of_elements) = tmp_chemical(1:fix_no_of_elements)
      fix_weighted_for_element(1:fix_no_of_elements) = tmp_chemical(1:fix_no_of_elements)
      if (weighted_3ch ) then 
        fix_weighted_for_element_3ch_ini(1: fix_no_of_elements ) = &
                                tmp_chemical(fix_no_of_elements+1: dim_chemical_function)
        fix_weighted_for_element_3ch(1: fix_no_of_elements ) = &
                                tmp_chemical(fix_no_of_elements+1: dim_chemical_function)
      end if   
    end if
    itopt = itopt + 1
    call main_train_mld
    !debug call status_allocate_config_desc
    mae_energy = train_mae_energy
    mae_force = train_mae_force
    mae_stress = train_mae_stress
    objval = factor_energy_error*mae_energy + factor_force_error*mae_force + factor_stress_error*mae_stress
    if (itopt == 1 ) then 
      if (rangml==0)   write(6,'("ML:                      itopt          obj               MAE_E              MAE_F               MAE_S")')
    end if   
    if (rangml==0) write (nunit, '(i7,4e20.10)') itopt, objval, mae_energy, mae_force, mae_stress
    if ((mod(itopt,10)==0).or.(itopt==1)) then 
    if (rangml==0) write (*,    '("  chem genetic --->  ", 1i7,4e20.10)') itopt,  objval, mae_energy, mae_force, mae_stress
    end if 
  end subroutine function_to_min_chemical 
  

end module hyper_optimize_chemical 


subroutine mld_optimize_weights()
  use ml_in_ndm_module, only: rangml
  use module_optimization, only:  max_iter_optimize_weights, optimize_ga_population, &
                                  dim_weights_function, dim_weights_function_full, &
                                  dim_chemical_function, & 
                                  optimize_weights_db, optimize_weights_chem, optimize_weights_Le, & 
                                  upper_weights, lower_weights, &
                                  tmp_chemical, upper_chemical, lower_chemical, itopt
  use data_type, only: nunit
  use mld_unit
  use mod_function, only: function_to_min
  use mld_logger
  use mld_string
  use mld_mpi
  use hyper_optimize_chemical, only: get_dim_objective_function_chemical, &
                                     init_from_chemical_in_objective_function, &
                                    function_to_min_chemical

  implicit none

  real(kind(0.d0))     :: VTR, F_XC, CR_XC, F_CR
  integer  :: NP, itermax, strategy, refresh, iwrite, nfeval
  !real(kind(0.d0)), dimension(dim_weights_function_full)  :: bestmem_XC
  real(kind(0.d0)), dimension(:), allocatable  :: bestmem_XC
  real(kind(0.d0))     :: bestval
  integer, dimension(3), parameter :: method = (/0, 1, 0/)
  integer  :: ic

  _NAMECURRENT_("mld_optimize_weights")

  _MLD_BEGIN_

  !external function_to_min

  ! size the dimension of J(w_1, w_2 ....)
  !debug if (rangml==0) write²(6,*) 'enter get_dim'
  if (optimize_weights_db) then 
    call get_dim_objective_function_weigths_database(dim_weights_function, dim_weights_function_full)

    call log_info('ML: genetic DB optimization for DB weights. dimensions partial '//vtoa(dim_weights_function)//nwl// &
    'and full  '//vtoa(dim_weights_function_full))

    ! From db put the values of w_1, w_2 ...
    call init_weigths_from_db_in_objective_function
    ! From w_1, w_2  fill the matrix of weights - the one, of dimension dim_data_train,
    ! Fill tmp_weights_snap the used by Amat in order to get the parameters
    if (allocated(bestmem_XC)) deallocate(bestmem_XC)
    allocate(bestmem_XC(dim_weights_function_full))
    itopt = 0 
  end if 

  if (optimize_weights_chem) then
    call get_dim_objective_function_chemical(dim_chemical_function)
    call log_info('ML: genetic DB optimization for chemical weights. dimension is  '//vtoa(dim_chemical_function))
    call init_from_chemical_in_objective_function(dim_chemical_function, tmp_chemical, upper_chemical, lower_chemical)
    if (allocated(bestmem_XC)) deallocate(bestmem_XC)
    allocate(bestmem_XC(dim_chemical_function))
    itopt = 0 
  end if 
  !the expected fitness value to reach.
  VTR = 1.d-04
  ! Population size.old was 40
  NP = optimize_ga_population
  ! The maximum number of iteration.
  itermax = max_iter_optimize_weights
  !Mutation scaling factor for real decision parameters.
  F_XC = 0.8d0
  !Crossover factor for real decision parameters.
  !was 0.8 
  CR_XC = 0.8d0
  !The strategy of the mutation operations is used in HDE 1 to 6
  strategy = 6
  !The unit specfier for writing to an external data file.
  iwrite = 7
  !The intermediate output will be produced after "refresh"
  !iterations. No intermediate output will be produced if
  !"refresh < 1".
  refresh = 200
  F_CR = 0.8d0
  !method(1) = 0, Fixed mutation scaling factors (F_XC)
  !          = 1, Random mutation scaling factors F_XC=[0, 1]
  !          = 2, Random mutation scaling factors F_XC=[-1, 1]
  !method(2) = 1, Random combined factor (F_CR) used for strategy = 6
  !               in the mutation operation
  !          = other, fixed combined factor provided by the user
  !method(3) = 1, Saving results in a data file.
  !          = other, displaying results only.
  open (newunit=nunit, file='optimization_error.dat', status='unknown')
  !debug if (rangml==0) write (6,*) 'Enter DE_Fortran'
  if (optimize_weights_db) then 
    call DE_Fortran90(function_to_min, dim_weights_function_full, lower_weights, upper_weights, &
                      VTR, NP, itermax, F_XC, CR_XC, strategy, refresh, iwrite, bestmem_XC, &
                      bestval, nfeval, F_CR, method)
    if (rangml == 0) write (6, '("ML: Number of function evaluation for genetic algo minimization ...:", i8)') nfeval
    if (rangml == 0) write (6, '("ML: The best value of the objective function ......................:", e20.10)') bestval
    if (rangml == 0) write (6, '("ML: The set of weigths ............................................:", i8)')
    do ic = 1, dim_weights_function
      if (rangml == 0) write (6, '(e20.10)') bestmem_XC(ic)
    end do
    if (rangml == 0) write (6, '("ML: The best lambda_krr L2 found .....................................:", e20.10)') bestmem_XC(dim_weights_function_full)
    if (optimize_weights_Le) then
      if (rangml == 0) write (6, '("ML: The best lambda_krr L1 found .....................................:", 2e20.10)') bestmem_XC(dim_weights_function_full - 1)
    else
    end if

    call ga_best_evaluation_db(bestmem_XC)

  end if 

  if (optimize_weights_chem) then 
    call DE_Fortran90(function_to_min_chemical, dim_chemical_function , lower_chemical, upper_chemical, &
                      VTR, NP, itermax, F_XC, CR_XC, strategy, refresh, iwrite, bestmem_XC, &
                      bestval, nfeval, F_CR, method)
    if (rangml == 0) write (6, '("ML: Number of function evaluation for genetic algo minimization ...:", i8)') nfeval
    if (rangml == 0) write (6, '("ML: The best value of the objective function ......................:", e20.10)') bestval
                                    
    call ga_best_evaluation_chemical(bestmem_XC)
  end if   


  _MLD_END_

end subroutine mld_optimize_weights




subroutine ga_best_evaluation_db(xval)
  use ml_in_ndm_module, only: lambda_krr
  use module_optimization, only: tmp_weights, dim_weights_function, dim_weights_function_full
  implicit none
  real(kind(0.d0)), dimension(dim_weights_function_full), intent(in)   :: xval

  tmp_weights(1:dim_weights_function) = xval(1:dim_weights_function)
  lambda_krr = xval(dim_weights_function_full)
  call fill_tmp_weigths_with_w
  call train_snap_get_parameters

end subroutine ga_best_evaluation_db

subroutine ga_best_evaluation_chemical(xval)
  use module_chemical_species, only: fix_ch_elements, fix_no_of_elements
  use ml_in_ndm_module, only: rangml, weighted_auto, weighted, weighted_3ch, & 
                              fix_weighted_for_element, fix_weighted_for_element_3ch, & 
                              fix_weighted_for_element_ini, fix_weighted_for_element_3ch_ini 
  use module_optimization, only: factor_energy_error, factor_force_error, factor_stress_error, & 
                                dim_chemical_function, tmp_chemical, itopt 
  use snap, only: train_mae_energy, train_mae_force, train_mae_stress  
  use module_test_mld, only: main_test_mld
  use main_mld_mod, only: main_train_mld


  use mld_logger
  implicit none
  real(kind(0.d0)), dimension(dim_chemical_function), intent(in)   :: xval
  integer :: ii, nunit 
  real(kind(0.d0))    :: objval
  real(kind(0.d0))     :: mae_energy, mae_force, mae_stress
  character(len=1)     :: quote, dquote
  quote = char(39)
  dquote = char(34)


  tmp_chemical(1:dim_chemical_function) = xval(1:dim_chemical_function)
  weighted_auto = .false. 
  if (weighted) then 
    fix_weighted_for_element_ini(1:fix_no_of_elements) = tmp_chemical(1:fix_no_of_elements)
    fix_weighted_for_element(1:fix_no_of_elements) = tmp_chemical(1:fix_no_of_elements)
    if (weighted_3ch ) then 
      fix_weighted_for_element_3ch_ini(1: fix_no_of_elements ) = &
                              tmp_chemical(fix_no_of_elements+1: dim_chemical_function)
      fix_weighted_for_element_3ch(1: fix_no_of_elements ) = &
                              tmp_chemical(fix_no_of_elements+1: dim_chemical_function)
    end if   
  end if

  if (rangml==0) then 
    open (file='optimization_chemical.info', newunit=nunit, action='write', status='unknown')
    write(nunit, '("chemical_elements = ")', advance='no') 
    write(nunit, '(a2)', advance="no") dquote 
    do ii = 1, fix_no_of_elements
      write(nunit, '(" ", a2)', advance="no") fix_ch_elements(ii)
    end do   
    write(nunit, '(" ", a2)') dquote  
   
    call log_info('ML: -----after genetic optimization the final weights ........')
    if (weighted) then
      call log_info('ML: 2nd channel ')

      do ii = 1, fix_no_of_elements
        write(*, '(a4, e20.7)') fix_ch_elements(ii), fix_weighted_for_element(ii)
      end do 
    
      write(nunit, '("weight_per_element = ")', advance="no") 
      write(nunit, '(a2)', advance="no") dquote 
      do ii =1, fix_no_of_elements
        write(nunit, '(" ", f10.4)', advance="no") fix_weighted_for_element(ii)
      end do   
      write(nunit, '(" ", a2)') dquote  

  
      if (weighted_3ch) then 
        call log_info('ML: 3rd channel ')
        do ii = 1, fix_no_of_elements
          write(*, '(a4, e20.7)') fix_ch_elements(ii), fix_weighted_for_element_3ch(ii)
        end do  

        write(nunit, '("weight_per_element_3ch = ")', advance="no") 
        write(nunit, '(a2)', advance="no") dquote 
        do ii =1, fix_no_of_elements
          write(nunit, '(" ", f10.4)', advance="no") fix_weighted_for_element_3ch(ii)
        end do   
        write(nunit, '(" ", a2)') dquote  
      end if    
    end if 
    close (nunit)
  end if ! rangml = 0 

  call main_train_mld
  mae_energy = train_mae_energy
  mae_force = train_mae_force
  mae_stress = train_mae_stress
  objval = factor_energy_error*mae_energy + factor_force_error*mae_force + factor_stress_error*mae_stress
  call log_info ('(ML: Final values for     itopt         obj               MAE_E              MAE_F               MAE_S :') 
  if (rangml==0) write (*,    '("  chem genetic --->  ", 1i7,4e20.10)') itopt,  objval, mae_energy, mae_force, mae_stress


end subroutine ga_best_evaluation_chemical





subroutine get_dim_objective_function_weigths_database(dim_weights_function, dim_weights_function_full)

  use ml_in_ndm_module, only: rangml
  use module_db_setup, only: iconf_data
  use snap, only:  map_weights_in_db, map_db_in_weights
  use module_optimization, only: optimize_weights_L1, optimize_weights_L2, optimize_weights_Le, &
                                 no_class_weights 
  use derived_types, only: db_model, config_real
  use mld_logger
  implicit none
  integer, intent(inout) :: dim_weights_function, dim_weights_function_full
  integer  :: i_w, idb, jc, i, iclass
  integer  :: ienergy, iforce, istress
  logical  :: l_energy, l_force, l_stress, skip_line

  l_energy = .false.
  l_force = .false.
  l_stress = .false.

  i_w = 0
  do idb = 1, size(db_model)
    !some classes are skipped
    skip_line = .false.
    do iclass = 1, size(no_class_weights)
      if (db_model(idb)%class == no_class_weights(iclass)) then
        skip_line = .true.
        db_model(idb)%has_optimize_weights_energy = .false.
        db_model(idb)%has_optimize_weights_force = .false.
        db_model(idb)%has_optimize_weights_stress = .false.
        cycle
      end if
    end do
    if (skip_line) cycle
    !end class selection

    ienergy = 0
    iforce = 0
    istress = 0
    !depending of fitting parameter T/F is in the db line the ienergy, iforce, istress is set to 1/0
    if (db_model(idb)%has_db_energy) then
      ienergy = 1
    end if
    if (db_model(idb)%has_db_force) then
      iforce = 1
    end if
    if (db_model(idb)%has_db_stress) then
      istress = 1
    end if

    !if all the componenets are false the ienergy, iforce. istress are putted to zero whatever is the initial value

    do jc = 1, iconf_data
      !selecting only trainning data ...
      if (.not. (config_real(jc)%train)) cycle
      ! in which line the jc conf is located ....
      if (config_real(jc)%db_line == idb) then
        l_energy = l_energy .or. config_real(jc)%has_energy
        l_force = l_force .or. config_real(jc)%has_force
        l_stress = l_stress .or. config_real(jc)%has_stress
      end if
    end do

    if (.not. l_energy) ienergy = 0
    if (.not. l_force) iforce = 0
    if (.not. l_stress) istress = 0
    i_w = i_w + ienergy + iforce + istress
    db_model(idb)%has_optimize_weights_energy = db_model(idb)%has_db_energy .and. l_energy
    db_model(idb)%has_optimize_weights_force = db_model(idb)%has_db_force .and. l_force
    db_model(idb)%has_optimize_weights_stress = db_model(idb)%has_db_stress .and. l_stress
    !if (rangml==0) write(*,*) idb, db_model(idb)%class, db_model(idb)%has_optimize_weights_energy, &
    !db_model(idb)%has_optimize_weights_force, db_model(idb)%has_optimize_weights_stress 
  end do


  dim_weights_function = i_w

  dim_weights_function_full = dim_weights_function
  if (optimize_weights_L1) dim_weights_function_full = dim_weights_function + 1
  if (optimize_weights_L2) dim_weights_function_full = dim_weights_function + 1
  if (optimize_weights_Le) dim_weights_function_full = dim_weights_function + 2

  if (allocated(map_weights_in_db)) deallocate (map_weights_in_db); allocate (map_weights_in_db(dim_weights_function))
  if (allocated(map_db_in_weights)) deallocate (map_db_in_weights); allocate (map_db_in_weights(size(db_model)))

  call log_info('ML: the dimension of the objective function, to optimize, is: '//vtoa(dim_weights_function))
  i_w = 0
  map_weights_in_db(:)%has_energy = .false. 
  map_weights_in_db(:)%has_force = .false. 
  map_weights_in_db(:)%has_stress = .false. 
  do i = 1, size(db_model)

    !some classes are skipped
    skip_line = .false.
    do iclass = 1, size(no_class_weights)
      if (db_model(i)%class == no_class_weights(iclass)) then
        skip_line = .true.
        cycle
      end if
    end do
    if (skip_line) cycle
    !end class selection


    if (db_model(i)%has_optimize_weights_energy) then
      i_w = i_w + 1
      map_weights_in_db(i_w)%db_line = i
      map_weights_in_db(i_w)%has_energy = .true.
      map_db_in_weights(i)%i_e = i_w
    end if
    if (db_model(i)%has_optimize_weights_force) then
      i_w = i_w + 1
      map_weights_in_db(i_w)%db_line = i
      map_weights_in_db(i_w)%has_force = .true.
      map_db_in_weights(i)%i_f = i_w
    end if
    if (db_model(i)%has_optimize_weights_stress) then
      i_w = i_w + 1
      map_weights_in_db(i_w)%db_line = i
      map_weights_in_db(i_w)%has_stress = .true.
      map_db_in_weights(i)%i_s = i_w
    end if
  end do

  if (i_w /= dim_weights_function) then
    if (rangml == 0) write (6, '("There are some problems in the definition of the dimension of objective function, old and new :", 2i6)') i_w, dim_weights_function
    stop 'dimension problem in get_dim_snap_objective_function'
  end if
end subroutine get_dim_objective_function_weigths_database




subroutine init_weigths_from_db_in_objective_function()

  use temporary_data_cov, only: dim_data_train
  use snap, only:  map_weights_in_db, tmp_weights_snap
  use derived_types, only: db_model
  use ml_in_ndm_module, only: rangml, lambda_krr, min_lambda_krr, max_lambda_krr
  use module_optimization, only: optimize_weights_Le, dim_weights_function, dim_weights_function_full, & 
                                 tmp_weights, upper_weights, lower_weights   
  use mld_logger
  implicit none
  integer  :: i
  real(kind(0.d0))     :: tmp_ini, tmp_fin

  if (allocated(tmp_weights_snap)) deallocate (tmp_weights_snap); allocate (tmp_weights_snap(dim_data_train))


  if (allocated(tmp_weights)) deallocate (tmp_weights); allocate (tmp_weights(dim_weights_function_full))
  if (allocated(lower_weights)) deallocate (lower_weights); allocate (lower_weights(dim_weights_function_full))
  if (allocated(upper_weights)) deallocate (upper_weights); allocate (upper_weights(dim_weights_function_full))
  do i = 1, dim_weights_function
    if (map_weights_in_db(i)%has_energy) then
      tmp_ini = db_model(map_weights_in_db(i)%db_line)%w_e_model
      tmp_fin = db_model(map_weights_in_db(i)%db_line)%w_e_end_model
      tmp_weights(i) = tmp_ini
      lower_weights(i) = tmp_ini
      upper_weights(i) = tmp_fin
      cycle 
    end if

    if (map_weights_in_db(i)%has_force) then
      tmp_ini = db_model(map_weights_in_db(i)%db_line)%w_f_model
      tmp_fin = db_model(map_weights_in_db(i)%db_line)%w_f_end_model
      tmp_weights(i) = tmp_ini
      lower_weights(i) = tmp_ini
      upper_weights(i) = tmp_fin
      cycle 
    end if


    if (map_weights_in_db(i)%has_stress) then
      tmp_ini = db_model(map_weights_in_db(i)%db_line)%w_s_model
      tmp_fin = db_model(map_weights_in_db(i)%db_line)%w_s_end_model
      tmp_weights(i) = tmp_ini
      lower_weights(i) = tmp_ini
      upper_weights(i) = tmp_fin
      cycle   
    end if


  end do



  ! tmp_weights(dim_weights_function+1)=(min_lambda_krr + max_lambda_krr)/2.d0
  if (optimize_weights_Le) then
    tmp_weights(dim_weights_function_full - 1) = min_lambda_krr
    lower_weights(dim_weights_function_full - 1) = min_lambda_krr
    upper_weights(dim_weights_function_full - 1) = max_lambda_krr


    tmp_weights(dim_weights_function_full) = min_lambda_krr
    lower_weights(dim_weights_function_full) = min_lambda_krr
    upper_weights(dim_weights_function_full) = max_lambda_krr

  else
    tmp_weights(dim_weights_function_full) = min_lambda_krr
    lower_weights(dim_weights_function_full) = min_lambda_krr
    upper_weights(dim_weights_function_full) = max_lambda_krr
  end if

  if (rangml == 0) then
  call log_info('ML:-------------------Genetic Optimization----------------------')
    if (lambda_krr > 0.d0) then
      write (6, *) 'ML: The optimization on weigths is performed in the same time with L2 norm min and max values if lambda: ', min_lambda_krr, max_lambda_krr
    end if
    write (6, *) 'ML: The size of the function to be optimized: ', dim_weights_function_full
  end if

end subroutine init_weigths_from_db_in_objective_function



subroutine fill_tmp_weigths_with_w
  use ml_in_ndm_module, only: rangml
  use module_db_setup, only: iconf_data
  use temporary_data_cov, only: dim_data_train
  use derived_types, only: config_real, db_model
  use snap, only: map_db_in_weights, tmp_weights_snap
  use module_optimization, only:  tmp_weights
  implicit none
  integer  :: i, itmp, itmp_f
  ! integer ::  jc
  ! logical :: lskip

  itmp = 0
  do i = 1, iconf_data
    if (.not. (config_real(i)%train)) cycle
    ! here is the problem

    !lskip=.false.
    !do jc=1,size(no_class_weights)
    !    if (config_real(i)%class==no_class_weights(jc)) then
    !         lskip=.true.
    !         cycle
    !       end if
    !end do
    !if (lskip) cycle

    if (config_real(i)%has_energy) then
      itmp = itmp + 1
      if (db_model(config_real(i)%db_line)%has_optimize_weights_energy) then
        tmp_weights_snap(itmp) = tmp_weights(map_db_in_weights(config_real(i)%db_line)%i_e)
      else
        tmp_weights_snap(itmp) = config_real(i)%w_e_model
      end if
    end if

    if (config_real(i)%has_force) then
      itmp_f = 3*config_real(i)%nat

      if (db_model(config_real(i)%db_line)%has_optimize_weights_force) then
        tmp_weights_snap(itmp + 1:itmp + itmp_f) = tmp_weights(map_db_in_weights(config_real(i)%db_line)%i_f)
      else
        tmp_weights_snap(itmp + 1:itmp + itmp_f) = config_real(i)%w_f_model
      end if
      itmp = itmp + itmp_f
    end if

    if (config_real(i)%has_stress) then
      if (db_model(config_real(i)%db_line)%has_optimize_weights_stress) then
        tmp_weights_snap(itmp + 1:itmp + 6) = tmp_weights(map_db_in_weights(config_real(i)%db_line)%i_s)
      else
        tmp_weights_snap(itmp + 1:itmp + 6) = config_real(i)%w_s_model
      end if
      itmp = itmp + 6
    end if
  end do

  if (itmp /= dim_data_train) then
    if (rangml == 0) write (6, '("There are some problems in tmp_weights_snap vector, old and new dimension:", 2i6)') dim_data_train, itmp
    stop 'dimension problem in fill_tmp_weigths_snap_with_w'
  end if
end subroutine fill_tmp_weigths_with_w



subroutine train_minimize_error(mae_energy, mae_force, mae_stress)

  use module_optimization, only: no_class_weights
  use snap, only: ene_snap, fit_snap                       
  use temporary_data_cov, only: yfunc_train, dim_data_train
  use derived_types, only: config_real
  use module_evaluate_parameters, only: product_w_params_Amat
  implicit none
  real(kind(0.d0)), intent(out)    :: mae_energy, mae_force, mae_stress
  integer  :: idata, i_w, jc, dim_y_force_error, dim_y_energy_error, dim_y_stress_error
  real(kind(0.d0))     :: force_snap, stress_snap
  !old things
  !character(len=2), dimension(2) :: class_energy_error
  !character(len=2), dimension(2) :: class_force_error
  real(kind(0.d0)), dimension(:), allocatable  :: y_energy_error_snap, y_energy_error_base, &
                                                  y_force_error_snap, y_force_error_base, &
                                                  y_stress_error_snap, y_stress_error_base
  real(kind(0.d0))     :: rmse_force, rmse_energy, rmse_stress
  !if (config_real(fit_snap(idata)%iconf)%class=="02") then
  !write (*,*) 'objective_test', idata,  config_real(fit_snap(idata)%iconf)%class, config_real(fit_snap(idata)%iconf)%filename
  !end if

  !class_energy_error(1)='02'
  !class_energy_error(2)='06'

  !class_force_error(1)='03'
  !class_force_error(2)='06'

  i_w = 0
  do idata = 1, dim_data_train
    if (fit_snap(idata)%energy) then

      do jc = 1, size(no_class_weights)
        if (config_real(fit_snap(idata)%iconf)%class == no_class_weights(jc)) cycle
        i_w = i_w + 1
      end do

    end if
  end do

  dim_y_energy_error = i_w


  i_w = 0
  do idata = 1, dim_data_train
    if (fit_snap(idata)%force) then

      do jc = 1, size(no_class_weights)
        if (config_real(fit_snap(idata)%iconf)%class == no_class_weights(jc)) cycle
        i_w = i_w + 1
      end do

    end if
  end do

  dim_y_force_error = i_w


  i_w = 0
  do idata = 1, dim_data_train
    if (fit_snap(idata)%stress) then

      do jc = 1, size(no_class_weights)
        if (config_real(fit_snap(idata)%iconf)%class == no_class_weights(jc)) cycle
        i_w = i_w + 1
      end do

    end if
  end do

  dim_y_stress_error = i_w



  if (allocated(y_force_error_base)) deallocate (y_force_error_base); allocate (y_force_error_base(dim_y_force_error))
  if (allocated(y_force_error_snap)) deallocate (y_force_error_snap); allocate (y_force_error_snap(dim_y_force_error))

  if (allocated(y_energy_error_base)) deallocate (y_energy_error_base); allocate (y_energy_error_base(dim_y_energy_error))
  if (allocated(y_energy_error_snap)) deallocate (y_energy_error_snap); allocate (y_energy_error_snap(dim_y_energy_error))


  if (allocated(y_stress_error_base)) deallocate (y_stress_error_base); allocate (y_stress_error_base(dim_y_stress_error))
  if (allocated(y_stress_error_snap)) deallocate (y_stress_error_snap); allocate (y_stress_error_snap(dim_y_stress_error))


  i_w = 0
  do idata = 1, dim_data_train
    if (fit_snap(idata)%energy) then
      do jc = 1, size(no_class_weights)
        if (config_real(fit_snap(idata)%iconf)%class == no_class_weights(jc)) cycle
        ! TODO ref energy missing ? 
        !ene_snap = dot_product(w_params(:, 1), Amat(:, idata))
        call product_w_params_Amat(ene_snap, idata)
        i_w = i_w + 1
        y_energy_error_snap(i_w) = ene_snap
        y_energy_error_base(i_w) = yfunc_train(idata)
      end do
    end if
  end do

  if (i_w > 0) then
    call rmse_mae(y_energy_error_base, y_energy_error_snap, dim_y_energy_error, rmse_energy, mae_energy)
  else
    rmse_force = 0.d0
    mae_force = 0.d0
  end if

  i_w = 0
  do idata = 1, dim_data_train
    if (fit_snap(idata)%force) then
      do jc = 1, size(no_class_weights)
        if (config_real(fit_snap(idata)%iconf)%class == no_class_weights(jc)) cycle
        !force_snap = dot_product(w_params(:, 1), Amat(:, idata))
        call product_w_params_Amat(force_snap, idata)
        i_w = i_w + 1
        y_force_error_snap(i_w) = force_snap
        y_force_error_base(i_w) = yfunc_train(idata)
      end do
    end if
  end do

  if (i_w > 0) then
    call rmse_mae(y_force_error_base, y_force_error_snap, dim_y_force_error, rmse_force, mae_force)
  else
    rmse_force = 0.d0
    mae_force = 0.d0
  end if


  i_w = 0
  do idata = 1, dim_data_train
    if (fit_snap(idata)%stress) then
      do jc = 1, size(no_class_weights)
        if (config_real(fit_snap(idata)%iconf)%class == no_class_weights(jc)) cycle
        !stress_snap = dot_product(w_params(:, 1), Amat(:, idata))
        call product_w_params_Amat(stress_snap, idata)
        i_w = i_w + 1
        y_stress_error_snap(i_w) = force_snap
        y_stress_error_base(i_w) = yfunc_train(idata)
      end do
    end if
  end do

  if (i_w > 0) then
    call rmse_mae(y_stress_error_base, y_stress_error_snap, dim_y_stress_error, rmse_stress, mae_stress)
  else
    rmse_force = 0.d0
    mae_force = 0.d0
  end if

end subroutine train_minimize_error




subroutine DE_Fortran90(obj, Dim_XC, XCmin, XCmax, VTR, NP, itermax, F_XC, &
                        CR_XC, strategy, refresh, iwrite, bestmem_XC, bestval, nfeval, &
                        F_CR, method)
  ! Berkley implementation from Taiwan
  !.......................................................................
  !
  ! Differential Evolution for Optimal Control Problems
  !
  !.......................................................................
  !  This Fortran 90 program translates from the original MATLAB
  !  version of differential evolution (DE). This FORTRAN 90 code
  !  has been tested on Compaq Visual Fortran v6.1.
  !  Any users new to the DE are encouraged to read the article of Storn and Price.
  !
  !  Refences:
  !  Storn, R., and Price, K.V., (1996). Minimizing the real function of the
  !    ICEC'96 contest by differential evolution. IEEE conf. on Evolutionary
  !    Comutation, 842-844.
  !
  !  This Fortran 90 program written by Dr. Feng-Sheng Wang
  !  Department of Chemical Engineering, National Chung Cheng University,
  !  Chia-Yi 621, Taiwan, e-mail: chmfsw@ccunix.ccu.edu.tw
  !.........................................................................
  !                obj : The user provided file for evlauting the objective function.
  !                      subroutine obj(xc,fitness)
  !                      where "xc" is the real decision parameter vector.(input)
  !                            "fitness" is the fitness value.(output)
  !             Dim_XC : Dimension of the real decision parameters.
  !      XCmin(Dim_XC) : The lower bound of the real decision parameters.
  !      XCmax(Dim_XC) : The upper bound of the real decision parameters.
  !                VTR : The expected fitness value to reach.
  !                 NP : Population size.
  !            itermax : The maximum number of iteration.
  !               F_XC : Mutation scaling factor for real decision parameters.
  !              CR_XC : Crossover factor for real decision parameters.
  !           strategy : The strategy of the mutation operations is used in HDE.
  !            refresh : The intermediate output will be produced after "refresh"
  !                      iterations. No intermediate output will be produced if
  !                      "refresh < 1".
  !             iwrite : The unit specfier for writing to an external data file.
  ! bestmen_XC(Dim_XC) : The best real decision parameters.
  !              bestval : The best objective function.
  !             nfeval : The number of function call.
  !         method(1) = 0, Fixed mutation scaling factors (F_XC)
  !                   = 1, Random mutation scaling factors F_XC=[0, 1]
  !                   = 2, Random mutation scaling factors F_XC=[-1, 1]
  !         method(2) = 1, Random combined factor (F_CR) used for strategy = 6
  !                        in the mutation operation
  !                   = other, fixed combined factor provided by the user
  !         method(3) = 1, Saving results in a data file.
  !                   = other, displaying results only.

  use data_type, only: IB, RP
  implicit none
  integer(kind=IB), intent(in)     :: NP, Dim_XC, itermax, strategy, &
                                      iwrite, refresh
  real(kind=RP), intent(in)  :: VTR, CR_XC
  real(kind=RP)  :: F_XC, F_CR
  real(kind=RP), dimension(Dim_XC), intent(in) :: XCmin, XCmax
  real(kind=RP), dimension(Dim_XC), intent(inout)    :: bestmem_XC
  real(kind=RP), intent(out) :: bestval
  integer(kind=IB), intent(out)    :: nfeval
  real(kind=RP), dimension(NP, Dim_XC)   :: pop_XC, bm_XC, mui_XC, mpo_XC, &
                                            popold_XC, rand_XC, ui_XC
  real(kind=RP), dimension(:), allocatable :: tmpvec
  integer(kind=IB)     :: i, ibest, iter
  integer(kind=IB), dimension(NP)  :: rot, a1, a2, a3, a4, a5, rt
  integer(kind=IB), dimension(4)   :: ind
  real(kind=RP)  :: tempval
  real(kind=RP), dimension(NP)     :: val
  real(kind=RP), dimension(Dim_XC) :: bestmemit_XC
  real(kind=RP), dimension(Dim_XC) :: rand_C1
  integer(kind=IB), dimension(3), intent(in)   :: method
  external obj
  intrinsic max, min, random_number, mod, abs, any, all, maxloc

  interface
    function randperm(num)
      use data_type, only: IB
      implicit none
      integer(kind=IB), intent(in)     :: num
      integer(kind=IB), dimension(num) :: randperm
    end function randperm
  end interface

  ! Initialize a population
  pop_XC = 0.0_RP
  do i = 1, NP
    call random_number(rand_C1)
    pop_XC(i, :) = XCmin + rand_C1*(XCmax - XCmin)
  end do


  ! Evaluate fitness functions and find the best member
  val = 0.0_RP
  nfeval = 0
  ibest = 1
  if (allocated(tmpvec)) deallocate(tmpvec) 
  allocate(tmpvec(size(pop_XC,2)))
  tmpvec(:) = pop_XC(1,:)
  call obj(tmpvec, val(1))
  !call obj(pop_XC(1, :), val(1))
  deallocate(tmpvec)
  bestval = val(1)
  nfeval = nfeval + 1
  do i = 2, NP
    if (allocated(tmpvec)) deallocate(tmpvec) 
    allocate(tmpvec(size(pop_XC,2)))
    tmpvec(:) = pop_XC(i,:)  
    !call obj(pop_XC(i, :), val(i))
    call obj(tmpvec, val(i))
    deallocate(tmpvec)
    nfeval = nfeval + 1
    if (val(i) < bestval) then
      ibest = i
      bestval = val(i)
    end if
  end do
  bestmemit_XC = pop_XC(ibest, :)
  bestmem_XC = bestmemit_XC

  bm_XC = 0.0_RP
  rot = (/(i, i=0, NP - 1)/)
  iter = 1

  ! Perform evolutionary computation
  do while (iter <= itermax)
    popold_XC = pop_XC

    ! Mutation operation
    ind = randperm(4)
    a1 = randperm(NP)
    rt = mod(rot + ind(1), NP)
    a2 = a1(rt + 1)
    rt = mod(rot + ind(2), NP)
    a3 = a2(rt + 1)
    rt = mod(rot + ind(3), NP)
    a4 = a3(rt + 1)
    rt = mod(rot + ind(4), NP)
    a5 = a4(rt + 1)
    bm_XC = spread(bestmemit_XC, DIM=1, NCOPIES=NP)

    ! Generating a random sacling factor
    select case (method(1))
    case (1)
      call random_number(F_XC)
    case (2)
      call random_number(F_XC)
      F_XC = 2.0_RP*F_XC - 1.0_RP
    end select

    ! select a mutation strategy
    select case (strategy)
    case (1)
      ui_XC = bm_XC + F_XC*(popold_XC(a1, :) - popold_XC(a2, :))

    case default
      ui_XC = popold_XC(a3, :) + F_XC*(popold_XC(a1, :) - popold_XC(a2, :))

    case (3)
      ui_XC = popold_XC + F_XC*(bm_XC - popold_XC + popold_XC(a1, :) - popold_XC(a2, :))

    case (4)
      ui_XC = bm_XC + F_XC*(popold_XC(a1, :) - popold_XC(a2, :) + popold_XC(a3, :) - popold_XC(a4, :))

    case (5)
      ui_XC = popold_XC(a5, :) + F_XC*(popold_XC(a1, :) - popold_XC(a2, :) + popold_XC(a3, :) &
                                       - popold_XC(a4, :))
    case (6)                ! A linear crossover combination of bm_XC and popold_XC
      if (method(2) == 1) call random_number(F_CR)
      ui_XC = popold_XC + F_CR*(bm_XC - popold_XC) + F_XC*(popold_XC(a1, :) - popold_XC(a2, :))

    end select

    ! Crossover operation
    call random_number(rand_XC)
    mui_XC = 0.0_RP
    mpo_XC = 0.0_RP
    where (rand_XC < CR_XC)
      mui_XC = 1.0_RP
      ! mpo_XC = 0.0_RP
    elsewhere
      ! mui_XC = 0.0_RP
      mpo_XC = 1.0_RP
    end where

    ui_XC = popold_XC*mpo_XC + ui_XC*mui_XC

    ! Evaluate fitness functions and find the best member
    do i = 1, NP
      ! Confine each of feasible individuals in the lower-upper bound
      ui_XC(i, :) = max(min(ui_XC(i, :), XCmax), XCmin)
      if (allocated(tmpvec)) deallocate(tmpvec) 
      allocate(tmpvec(size(ui_XC,2)))
      tmpvec(:) = ui_XC(i,:)  
      call obj(tmpvec, tempval)
      deallocate(tmpvec)
      !call obj(ui_XC(i, :), tempval)
      nfeval = nfeval + 1
      if (tempval < val(i)) then
        pop_XC(i, :) = ui_XC(i, :)
        val(i) = tempval
        if (tempval < bestval) then
          bestval = tempval
          bestmem_XC = ui_XC(i, :)
        end if
      end if
    end do
    bestmemit_XC = bestmem_XC
    if ((refresh > 0) .and. (mod(iter, refresh) == 0)) then
      if (method(3) == 1) write (unit=iwrite, FMT=203) iter
      write (unit=*, FMT=203) iter
      do i = 1, Dim_XC
        if (method(3) == 1) write (unit=iwrite, FMT=202) i, bestmem_XC(i)
        write (*, FMT=202) i, bestmem_XC(i)
      end do
      if (method(3) == 1) write (unit=iwrite, FMT=201) bestval
      write (unit=*, FMT=201) bestval
    end if
    iter = iter + 1
    if (bestval <= VTR .and. refresh > 0) then
      write (unit=iwrite, FMT=*) ' The best fitness is smaller than VTR'
      write (unit=*, FMT=*) 'The best fitness is smaller than VTR'
      exit
    end if
  end do
  ! end the evolutionary computation

201 format(2x, 'bestval =', ES14.7,/)
202 format(5x, 'bestmem_XC(', I3, ') =', ES12.5)
203 format(2x, 'No. of iteration  =', I8)

end subroutine DE_Fortran90



