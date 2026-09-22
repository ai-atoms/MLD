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



subroutine md_init_potential_ml()
  ! MD drivers for MLD! Should be updated. 
  ! Here the ml potential is initialized. This subroutine is used
  ! by MD program and is called in NDM's init.F90

#ifdef MLD_NDM
  use gen_com_m, only: umass, dmtype
  use gen_com_m_ml, only: rvois
  use var_pot, only: cm, ntyp, typ_and_pot, typ_pot_pair, lpotentiel, rue_pot, rue_pair, &
                           ipotentiel, npotmax, iewald, l3c, r3cm, rumax, npair
#else
  use ondm_gen_com_m, only: rvois, umass, dmtype, A2cm
  use ondm_var_pot, only: cm, ntyp, typ_and_pot, typ_pot_pair, lpotentiel, rue_pot, rue_pair, &
                           ipotentiel, npotmax, iewald, l3c, r3cm, rumax, npair
#endif
  use time_measure, only: temps_energy, temps_force, temps_descripteurs, temps_neigh, temps_stress, init_time_measure
  use ml_in_ndm_module, only: debug, prepare_factorial, toy_model, ML_MLD_DMTYPE, &
                    ml_type, ml_type_krr
  use module_neigh_local, only: r_cut 
  use module_chemical_species, only:  periodic_table_element, fix_type_to_periodic, fix_no_of_elements
  use module_ml_scalapack, only: scalapack_driver
  use main_mld_mod, only: read_parameters_for_md, set_mld_dimension
  use module_ftnbody, only: ftnbody_read_from_potential
  use module_kernel_2b, only: activate_k2b, dim_kernel_2b
  use module_kernel, only: kernel_type, kernel_random, kernel_random_po, kernel_random_maha
  use snap, only: w_params, dim_design_line
  use mld_logger
  use mld_string
  use mld_mpi

  implicit none

  _NAMECURRENT_("md_init_potential_ml")


  integer  :: i


  _MLD_BEGIN_
  if (debug) call log_info("... enter in : "//NAMECURRENT)

  !HERE ntyp should be changed ....

  call read_ml_file()

  call log_info("ML: read_ml_file in MD_mode ... done ")
  ntyp = fix_no_of_elements
  call alloc_typ_ml

  call log_info("ML: alloc_typ_ml in MD_mode ... done ")

  !HERE ntyp and cm set by hand ...
  !MAYBE cm should be initialised from ml file ...
  !TODOmd put the chemical_symbols_somewhere.
  do i = 1, ntyp
    call log_warning("ML: MD_mode Chemical and type are not fixed in MD mode. More than one type nyp = 1 is DANGEROUS")
    cm(:ntyp) = periodic_table_element(fix_type_to_periodic(i))%mass
  end do

  if ((dmtype /= ML_MLD_DMTYPE) .and. (ntyp >= 2)) then
    call log_warning('ML: MD_mode other WARNING  for ntyp > 1 not yet implementation for many elements systems')
    call log_warning('ML: MD_mode TODOmd')
  end if

  !Fe mass
  !cm(:ntyp)=periodic_table_element(26)%mass
  cm(:ntyp) = cm(:ntyp)*umass
  call log_info('ML: MD_mode initial type and masses in uam '//vtoa(ntyp)//nwl// 'mass setting '//vtoa(cm(:ntyp)/umass))

  !here im not sure. Maybe rvois should be read in * .din and rue_pot fixed to r_cut of desciptors
  !like that rue_pot(:)=r_cut*A2cm

  rvois=r_cut*A2cm
  rue_pot(:) = rvois

  do i = 1, npair
    if (typ_pot_pair(i) == ipotentiel) rue_pair(i) = rue_pot(ipotentiel)
  end do

  iewald = 0; l3c = .false.; r3cm = 0.d0
  allocate (typ_and_pot(ntyp, npotmax))
  typ_and_pot(1:ntyp, ipotentiel) = .true.
  typ_pot_pair(:) = ipotentiel

  if (ipotentiel /= 20) call mld_mpi_abort('ML: MD_mode ipotentiel should be 20 in din file')

  lpotentiel(ipotentiel) = .true.
  rue_pot(ipotentiel) = rvois
  rue_pair(:) = rue_pot(ipotentiel)
  rumax = rvois
  !initialisation ... very general

  if (dmtype /= ML_MLD_DMTYPE) then
    call prepare_factorial()

    call periodic_table()
    call log_info("ML: MD_mode periodic_table ... done")
    if (scalapack_driver) then 
      call set_up_scalapack
      call log_info("ML: MD_mode  set_up_scalapack  ... done")
    end if

    if (toy_model) call set_toy_model
    ! MD / prediction: load random-sampling descriptors (ftnbody) from random.xml
    ! instead of re-estimating them; init_descriptors honours this flag.
    ftnbody_read_from_potential = .true.
    call init_descriptors
    call log_info("ML: MD_mode  init_descriptors   ... done")
    call log_warning("TODOmd There are k2b, kernels ... which are not initialized in this MD brach.")
     ! init the kernel if there is some ...
     dim_kernel_2b = 0
     if (activate_k2b) then 
       call init_kernel_k2b
       call log_info('ML: MD_mode 2-body kernel is activated and the dimension is '//vtoa(dim_kernel_2b))
     end if 

    if (ml_type == ml_type_krr) then
      if (kernel_type == kernel_random) then
        call init_sample_kernel
      else if (kernel_type == kernel_random_po) then
        call init_sample_kernel_po
      else if (kernel_type == kernel_random_maha) then 
        call init_sample_kernel_maha  
      else
        call init_kernel
      end if
    end if


    call log_warning("TODOmd k2, kernels etc things to do in md_init_potential_ml")
    !initialisation ... only for snap
    !call md_allocate_snap_params()
    call set_mld_dimension
    call log_info("ML: MD_mode  set_mld_dimension   ... done")

    if (allocated(w_params)) deallocate (w_params); allocate (w_params(dim_design_line, 1))
    w_params(:, 1) = 0.d0
    if (mld_rank==0) call read_parameters_for_md
    call comm_mld%bcast(0, w_params(:, 1))
    call log_info("ML: MD_mode read_parameters_for_md ... done") 
  end if 

  !initialization of time counter ...
  call init_time_measure()
  temps_force = 0.d0
  temps_neigh = 0.d0
  temps_energy = 0.d0
  temps_descripteurs = 0.d0
  temps_stress = 0.d0
  if (debug) call log_info("... exit from: "//NAMECURRENT)

  _MLD_END_
end subroutine md_init_potential_ml



!NoMoreInUse! subroutine md_init_config_ml()
!NoMoreInUse!   ! MD_mode drivers for MLD! Should be updated. 
!NoMoreInUse!   ! Here the ml potential is initialized. This subroutine is used
!NoMoreInUse!   ! by MD program and is called in NDM's init.F90 sub after that configuration was read
!NoMoreInUse! 
!NoMoreInUse!   use ml_in_ndm_module, only: write_desc, ML_MLD_DMTYPE
!NoMoreInUse!   use module_db_setup, only: md_iconf, iconf_data, iconf_data_train, iconf_data_test
!NoMoreInUse! #ifdef MLD_NDM
!NoMoreInUse!   use gen_com_m, only: A2cm, dmtype, fnam
!NoMoreInUse!   use gen_com_m_ml, only: im, imm
!NoMoreInUse!   !use var_pot, only: ntyp
!NoMoreInUse!   use ondm_gen_com_m, only: volu
!NoMoreInUse! #else
!NoMoreInUse!   use ondm_gen_com_m, only: im, imm,  dmtype, fnam
!NoMoreInUse!   !use ondm_var_pot, only: ntyp
!NoMoreInUse! #endif
!NoMoreInUse!   use derived_types, only: config_real, config_desc 
!NoMoreInUse!   !use module_db_poscar, only: fix_ref_energy_per_element
!NoMoreInUse!   use module_md_mld, only: md_allocate_mld_desc
!NoMoreInUse!   use snap, only: fp_snap
!NoMoreInUse! 
!NoMoreInUse!   use main_mld_mod
!NoMoreInUse!   use mld_logger
!NoMoreInUse!   use mld_string
!NoMoreInUse!   use mld_mpi
!NoMoreInUse! 
!NoMoreInUse!   implicit none
!NoMoreInUse!   !integer  :: nb_elements
!NoMoreInUse!   character(len=:), allocatable    :: fnamposcar 
!NoMoreInUse!   logical :: poscarok 
!NoMoreInUse! 
!NoMoreInUse!   _NAMECURRENT_("md_init_config_ml")
!NoMoreInUse! 
!NoMoreInUse!   _MLD_BEGIN_
!NoMoreInUse!   call log_warning("TODOmd ....") 
!NoMoreInUse!   call log_warning("md_ini_config_ml: Here should be implemented read of some poscar configuration")
!NoMoreInUse!   call log_warning("md_ini_config_ml: For the moment is just gin reading. From disk or HARDGIN  mode")
!NoMoreInUse! 
!NoMoreInUse!   !initialisation ... very general
!NoMoreInUse!   iconf_data = 1
!NoMoreInUse!   iconf_data_train = 1
!NoMoreInUse!   iconf_data_test = 0
!NoMoreInUse!   ! config_real and condif_desc objects
!NoMoreInUse!   if (allocated(config_real)) deallocate (config_real); allocate (config_real(iconf_data))
!NoMoreInUse!   if (allocated(config_desc)) deallocate (config_desc); allocate (config_desc(iconf_data))
!NoMoreInUse! 
!NoMoreInUse!   md_iconf = iconf_data
!NoMoreInUse!   if (imm .gt. im) then
!NoMoreInUse!     call log_warning('imm should be resized to im. MILADY cannot work otherwise '//vtoa([im, imm]))
!NoMoreInUse!   end if
!NoMoreInUse! 
!NoMoreInUse!   call md_allocate_mld_desc(fp_snap)
!NoMoreInUse! 
!NoMoreInUse! 
!NoMoreInUse!   fnamposcar = fnam//'.poscar'
!NoMoreInUse!   inquire (file=fnamposcar, exist=poscarok)
!NoMoreInUse!   if (poscarok) then 
!NoMoreInUse!      config_real(md_iconf)%filename = fnamposcar 
!NoMoreInUse!      call log_info('ML:  MD_mode poscar file found: '//fnamposcar)
!NoMoreInUse!   else 
!NoMoreInUse!      call log_critical('ML:  MD_mode poscar file not found. Please provide a poscar file as input')
!NoMoreInUse!      call log_critical('ML:  We exxpect a file with the name: '//fnamposcar)
!NoMoreInUse!      call mld_critical_abort('ML:  MD_mode ...  please provide a poscar file as input')
!NoMoreInUse!   end if    
!NoMoreInUse! 
!NoMoreInUse! 
!NoMoreInUse!   !config_real(md_iconf)%volume = volu/A2cm**3
!NoMoreInUse!   config_real(md_iconf)%has_energy = .true.
!NoMoreInUse!   config_real(md_iconf)%has_force = .true.
!NoMoreInUse!   config_real(md_iconf)%has_stress = .true.
!NoMoreInUse! 
!NoMoreInUse! 
!NoMoreInUse!   if (write_desc .and. (dmtype /= ML_MLD_DMTYPE)) then
!NoMoreInUse!     call log_warning('ML: writing descriptors using MD is not possible yet. This option will be switch to false')
!NoMoreInUse!     write_desc = .false.
!NoMoreInUse!   end if
!NoMoreInUse! 
!NoMoreInUse!   !$! !NOLD ... in all next lines was im instead of imm ...
!NoMoreInUse!   !$! if (allocated(config_real(md_iconf)%itype)) deallocate (config_real(md_iconf)%itype)
!NoMoreInUse!   !$! allocate (config_real(md_iconf)%itype(imm))
!NoMoreInUse!   !$! if (allocated(config_real(md_iconf)%pos_cart)) deallocate (config_real(md_iconf)%pos_cart)
!NoMoreInUse!   !$! allocate (config_real(md_iconf)%pos_cart(3, imm))
!NoMoreInUse!   !$! if (allocated(config_real(md_iconf)%pos_crst)) deallocate (config_real(md_iconf)%pos_crst)
!NoMoreInUse!   !$! allocate (config_real(md_iconf)%pos_crst(3, imm))
!NoMoreInUse!   !$! if (allocated(config_real(md_iconf)%force)) deallocate (config_real(md_iconf)%force)
!NoMoreInUse!   !$! allocate (config_real(md_iconf)%force(3, imm))
!NoMoreInUse!   !$! if (allocated(config_real(md_iconf)%mass_per_type)) deallocate (config_real(md_iconf)%mass_per_type)
!NoMoreInUse!   !$! allocate (config_real(md_iconf)%mass_per_type(ntyp))  
!NoMoreInUse!   !$! if ((dmtype /= ML_MLD_DMTYPE) .and. (ntyp >= 2)) then
!NoMoreInUse!   !$!   call log_warning('not yet implementation for many elements systems (dmtype, ntyp)')
!NoMoreInUse!   !$! end if  
!NoMoreInUse!   !$! nb_elements = ntyp
!NoMoreInUse!   !$! if (allocated(config_real(md_iconf)%weight_per_type)) deallocate (config_real(md_iconf)%weight_per_type)
!NoMoreInUse!   !$! allocate (config_real(md_iconf)%weight_per_type(nb_elements))
!NoMoreInUse!   !$! if (allocated(config_real(md_iconf)%weight_per_type_3ch)) deallocate (config_real(md_iconf)%weight_per_type_3ch)
!NoMoreInUse!   !$! allocate (config_real(md_iconf)%weight_per_type_3ch(nb_elements))
!NoMoreInUse!   !$! if (allocated(config_real(md_iconf)%ref_energy_per_element)) deallocate (config_real(md_iconf)%ref_energy_per_element)
!NoMoreInUse!   !$! allocate (config_real(md_iconf)%ref_energy_per_element(nb_elements))
!NoMoreInUse!   !$! if (allocated(config_real(md_iconf)%mass_per_type)) deallocate (config_real(md_iconf)%mass_per_type)
!NoMoreInUse!   !$! allocate (config_real(md_iconf)%mass_per_type(nb_elements))
!NoMoreInUse!   !$! if (allocated(config_real(md_iconf)%Z_per_type)) deallocate (config_real(md_iconf)%Z_per_type)
!NoMoreInUse!   !$! allocate (config_real(md_iconf)%Z_per_type(nb_elements))
!NoMoreInUse!   !$! if (allocated(config_real(md_iconf)%covalent_radius_per_type)) deallocate (config_real(md_iconf)%covalent_radius_per_type)
!NoMoreInUse!   !$! allocate (config_real(md_iconf)%covalent_radius_per_type(nb_elements))
!NoMoreInUse!   !$! if (allocated(config_real(md_iconf)%fix_type_poscar_to_periodic)) deallocate (config_real(md_iconf)%fix_type_poscar_to_periodic)
!NoMoreInUse!   !$! allocate (config_real(md_iconf)%fix_type_poscar_to_periodic(nb_elements))
!NoMoreInUse! 
!NoMoreInUse!   !$! TODOmd   MD_mode  type for MD - properly nowQD
!NoMoreInUse!   !$! do ij = 1, nb_elements
!NoMoreInUse!   !$!   icnt = 0
!NoMoreInUse!   !$!   do i_p = 1, fix_no_of_elements
!NoMoreInUse!   !$!     icnt = icnt + 1
!NoMoreInUse!   !$!     config_real(md_iconf)%ref_energy_per_element(ij) = fix_ref_energy_per_element(i_p)
!NoMoreInUse!   !$!   end do
!NoMoreInUse!   !$! end do
!NoMoreInUse! 
!NoMoreInUse! 
!NoMoreInUse!   _MLD_END_
!NoMoreInUse! end subroutine md_init_config_ml



subroutine put_ndm_into_ml_config(iconf)
   ! MD drivers for MLD! Should be updated. 

#ifdef MLD_NDM  
  use gen_com_m, only: A2cm, erg2ev, umass
  use gen_com_m_ml, only: im, imm, at, bg
  use ondm_gen_com_m, only: volu
  use tab_imm_m_ml, only: xp, fp, ityp
  use var_pot, only: ntyp, cm
#else
  use ondm_gen_com_m, only: im, at, bg, imm, A2cm, erg2ev, umass, volu
  use ondm_tab_imm_m, only: xp, fp, ityp
  use ondm_var_pot, only: ntyp, cm
#endif

  use derived_types, only: config_real
  use mld_logger

  implicit none

  _NAMECURRENT_("put_ndm_into_ml_config")


  integer, intent(in)  :: iconf

  _MLD_BEGIN_
  config_real(iconf)%volume = volu/A2cm**3
  config_real(iconf)%ntypes = ntyp
  config_real(iconf)%nat = im
  config_real(iconf)%im = im
  config_real(iconf)%imm = imm
  !write (6,*) 'im, imm', im, imm, size(ityp), size(config_real(iconf)%itype)
  !NOLD im=imm
  config_real(iconf)%itype(1:imm) = ityp(1:imm)
  config_real(iconf)%pos_cart(1:3, 1:imm) = xp(1:3, 1:imm)/A2cm
  !config_real(iconf)%pos_crst(3,1:imm) = xc(3,1:imm)/A2cm
  config_real(iconf)%force(1:3, 1:imm) = fp(1:3, 1:imm)*(A2cm*erg2ev)
  config_real(iconf)%cell = at(:, :)/A2cm
  config_real(iconf)%bg_cell = bg(:, :)*A2cm
  config_real(iconf)%mass_per_type(:) = cm(:)/umass
  !future config_real(iconf)%prev_pos_cart(1:3,1:imm)=xpp(1:3,1:imm)/A2cm
  _MLD_END_
end subroutine put_ndm_into_ml_config


!$! subroutine put_ml_config_into_ndm(iconf)
!$! 
!$!   use ml_in_ndm_module, only: debug, rangml
!$! #ifdef MLD_NDM
!$!   use gen_com_m, only: A2cm, erg2ev, umass
!$!   use gen_com_m_ml, only: im, imm, at, bg
!$!   use ondm_gen_com_m, only: volu
!$!   use var_pot, only: ntyp, cm
!$!   use tab_imm_m_ml, only: xp, fp, ityp
!$! #else
!$!   use ondm_gen_com_m, only: at, bg, im, imm, A2cm, erg2ev, umass, volu
!$!   use ondm_var_pot, only: ntyp, cm
!$!   use ondm_tab_imm_m, only: xp, fp, ityp
!$! #endif
!$!   use derived_types, only: config_real
!$!   use mld_logger
!$! 
!$!   implicit none
!$! 
!$!   _NAMECURRENT_("put_ml_config_into_ndmu")
!$! 
!$! 
!$!   integer, intent(in)  :: iconf
!$! 
!$! 
!$!   _MLD_BEGIN_
!$!   !NOLD imm = config_real(iconf)%nat
!$!   !THOSE TWO LINES ARE NOLD
!$!   im = config_real(iconf)%nat
!$!   imm = config_real(iconf)%imm
!$!   !NOLD END............
!$!   ityp(1:imm) = config_real(iconf)%itype(1:imm)
!$!   xp(1:3, 1:imm) = A2cm*config_real(iconf)%pos_cart(1:3, 1:imm)
!$!   !xc(3,1:imm)=A2cm*config_real(iconf)%pos_crst(3,1:imm)
!$!   fp(1:3, 1:imm) = config_real(iconf)%force(1:3, 1:imm)/(A2cm*erg2ev)
!$!   at(:, :) = config_real(iconf)%cell(:, :)*A2cm
!$!   bg(:, :) = config_real(iconf)%bg_cell(:, :)/A2cm
!$!   cm(:) = config_real(iconf)%mass_per_type(:)*umass
!$!   volu = config_real(iconf)%volume*A2cm**3
!$!   ntyp = config_real(iconf)%ntypes
!$!   !future xpp(1:3,1:imm)=config_real(iconf)%prev_pos_cart(1:3,1:imm)*A2cm
!$!   _MLD_END_
!$! end subroutine put_ml_config_into_ndm


subroutine md_calfo_ml
  ! MD drivers for MLD! Should be updated. 
  ! provides the energy and forces for NDM MD main program.
  ! called in the main NDM's calfo.F90
  ! Output: (to be completed)
  !         potist and fp ? through NDM module


  !use mpi
  !use mod_mpi_ml

#ifdef MLD_NDM
  use gen_com_m, only: potist, sig, ev2erg, erg2ev, A2cm, evA2dyn
  use tab_imm_m_ml, only: fp
#else
  use ondm_gen_com_m, only: potist, sig, ev2erg, erg2ev, A2cm, evA2dyn
  use ondm_tab_imm_m, only: fp
#endif

  use ml_in_ndm_module, only: ml_type, ml_type_basis, ml_type_krr, &
                              prepare_factorial
  use module_db_setup, only: md_iconf                            
  use snap, only: ene_snap, fp_snap, stress_snap
  use main_mld_mod
  use time_measure, only: temps_energy, temps_force, temps_descripteurs, temps_neigh, temps_stress
  use mld_logger
  use mld_mpi
  use mld_unit

  implicit none

  real(kind(0.d0))     :: temps1, temps2, temps3
  logical  :: post_desc
    _NAMECURRENT_("md_calfo_ml")


  _MLD_BEGIN_
  if ((ml_type == ml_type_basis) .or. (ml_type == ml_type_krr)) then
    !ph  if (rangml==0) write (*,*) 'here0', md_iconf, nox, noxyz
    call put_ndm_into_ml_config(md_iconf)
    !if (debug) write (6,*) 'was put_ml_config_into_ndm'
    !ph  if (rangml==0) write (*,*) 'here1', nox, noxyz
    call test_if_config_is_small(md_iconf)
    !  config_real(md_iconf)%small=.true.
    !ph  if (rangml==0) write (*,*) 'here2', nox,noxyz
    !if (debug) write (6,*) 'was test_if_config_is_small'
    temps1 = MPI_Wtime()
    call calc_neighbours(md_iconf)
    temps2 = MPI_Wtime()

    temps_neigh = temps_neigh + (temps2 - temps1)

    temps1 = MPI_Wtime()
    post_desc=.true. 
    call compute_descriptors(md_iconf, post_desc)

    temps2 = MPI_Wtime()
    temps_descripteurs = temps_descripteurs + (temps2 - temps1)

    temps1 = MPI_Wtime()

    call md_mld_compute_energy(md_iconf)

    temps2 = MPI_Wtime()
    call md_mld_compute_force(md_iconf)
    temps3 = MPI_Wtime()
    temps_energy = temps_energy + (temps2 - temps1)
    temps_force = temps_force + (temps3 - temps2)


    temps1 = MPI_Wtime()
    call md_mld_compute_stress(md_iconf)
    temps2 = MPI_Wtime()

    temps_stress = temps_stress + (temps2 - temps1)
    !put energy in NDM units
    potist = ene_snap*ev2erg
    !put energy in NDM units
    fp(:, :) = fp_snap(:, :)/(A2cm*erg2ev)
    sig(1, 1) = stress_snap(1)
    sig(2, 2) = stress_snap(2)
    sig(3, 3) = stress_snap(3)
    sig(2, 3) = stress_snap(4)
    sig(1, 3) = stress_snap(5)
    sig(1, 2) = stress_snap(6)
    sig(2, 1) = sig(1, 2)
    sig(3, 1) = sig(1, 3)
    sig(3, 2) = sig(2, 3)

    !put stres in NDM units
    sig(:, :) = sig(:, :)*1.d+09/evA2dyn
    !write (*,*) 'ENE_SNAP', ene_snap, potist, maxval(fp_snap)
    !stop 'ENE_SNAP00'
    !write (*,*) 'here5',  nox, noxyz, rumax
    !stop "testing clafo_ml"

  end if                  ! ml_type==ml_type_basis OR ml_type == ml_type_krr

  !dmtype_fake if (dmtype == 5) then
  !dmtype_fake   if (rangml == 0) then
  !dmtype_fake 
  !dmtype_fake     open (newunit=unitco, file='COORD', status='unknown')
  !dmtype_fake     open (newunit=unitfo, file='FORCE', status='unknown')
  !dmtype_fake 
  !dmtype_fake     ! Writting COORD
  !dmtype_fake     write (unitco, '(a)') 'ITEM: TIMESTEP'
  !dmtype_fake     write (unitco, '(a)') '0'
  !dmtype_fake     write (unitco, '(a)') 'ITEM: NUMBER OF ATOMS'
  !dmtype_fake     write (unitco, '(i6)') config_real(md_iconf)%nat
  !dmtype_fake     write (unitco, '(a)') 'ITEM: BOX BOUNDS xy xz yz pp pp pp'
  !dmtype_fake     write (unitco, '(3e23.15)') 0.d0, config_real(md_iconf)%cell(1, 1), 0.d0
  !dmtype_fake     write (unitco, '(3e23.15)') 0.d0, config_real(md_iconf)%cell(2, 2), 0.d0
  !dmtype_fake     write (unitco, '(3e23.15)') 0.d0, config_real(md_iconf)%cell(3, 3), 0.d0
  !dmtype_fake     tmp_val = config_real(md_iconf)%cell(1, 2)**2 + config_real(md_iconf)%cell(1, 3)**2 + &
  !dmtype_fake               config_real(md_iconf)%cell(2, 3)**2 + config_real(md_iconf)%cell(2, 1)**2 + &
  !dmtype_fake               config_real(md_iconf)%cell(3, 1)**2 + config_real(md_iconf)%cell(3, 2)**2
  !dmtype_fake     if (tmp_val >= 1.d-15) then
  !dmtype_fake       call mld_mpi_abort('this save type is not implemented for triclinic box. Stop in calfo_ml')
  !dmtype_fake     end if
  !dmtype_fake     write (unitco, '(a)') 'ITEM: ATOMS id xu yu zu'
  !dmtype_fake     do ia = 1, config_real(md_iconf)%nat
  !dmtype_fake       write (unitco, '(i5,3f30.15)') ia, config_real(md_iconf)%pos_cart(:, ia)
  !dmtype_fake     end do
  !dmtype_fake     close (unitco, status='keep')
  !dmtype_fake 
  !dmtype_fake     ! Writing FORCE
  !dmtype_fake   end if
  !dmtype_fake 
  !dmtype_fake   ! Writting COORD
  !dmtype_fake   write (unitfo, '(a)') 'ITEM: TIMESTEP'
  !dmtype_fake   write (unitfo, '(a)') '0'
  !dmtype_fake   write (unitfo, '(a)') 'ITEM: NUMBER OF ATOMS'
  !dmtype_fake   write (unitfo, '(i6)') config_real(md_iconf)%nat
  !dmtype_fake   write (unitfo, '(a)') 'ITEM: BOX BOUNDS xy xz yz pp pp pp'
  !dmtype_fake   write (unitfo, '(3e23.15)') 0.d0, config_real(md_iconf)%cell(1, 1), 0.d0
  !dmtype_fake   write (unitfo, '(3e23.15)') 0.d0, config_real(md_iconf)%cell(2, 2), 0.d0
  !dmtype_fake   write (unitfo, '(3e23.15)') 0.d0, config_real(md_iconf)%cell(3, 3), 0.d0
  !dmtype_fake   tmp_val = config_real(md_iconf)%cell(1, 2)**2 + config_real(md_iconf)%cell(1, 3)**2 + &
  !dmtype_fake             config_real(md_iconf)%cell(2, 3)**2 + config_real(md_iconf)%cell(2, 1)**2 + &
  !dmtype_fake             config_real(md_iconf)%cell(3, 1)**2 + config_real(md_iconf)%cell(3, 2)**2
  !dmtype_fake   if (tmp_val >= 1.d-15) then
  !dmtype_fake     call mld_mpi_abort('this save type is not implemented for triclinic box. Stop in calfo_ml')
  !dmtype_fake   end if
  !dmtype_fake   write (unitfo, '(a)') 'ITEM: ATOMS id fx fy fz '
  !dmtype_fake   do ia = 1, config_real(md_iconf)%nat
  !dmtype_fake     write (unitfo, '(i5,3f30.15)') ia, fp_snap(:, ia)
  !dmtype_fake   end do
  !dmtype_fake   close (unitfo, status='keep')
  !dmtype_fake   call mld_mpi_abort('TEST FORCES IN ML')
  !dmtype_fake end if

  _MLD_END_

end subroutine md_calfo_ml


!remove! subroutine md_allocate_snap_params()
!remove!    ! MD drivers for MLD! Should be updated. 
!remove!    ! used olny for md
!remove!    use ml_in_ndm_module, only:  mld_order, mld_linear, mld_quadratic
!remove!    use temporary_data_cov, only: dim_xdesc
!remove!    use snap, only: w_params
!remove!    use module_mld_quadratic, only: dim_xdesc_quadratic
!remove!    use mld_logger
!remove!    implicit none
!remove!    _NAMECURRENT_("md_allocate_snap_params")
!remove!    _MLD_BEGIN_
!remove!    if (dim_xdesc == 0) then
!remove!       call log_critical("ML: dim_xdesc Fatal Error")
!remove!       stop "dim_xdesc is zero in md_allocate_snap"
!remove!    end if
!remove! 
!remove!   call log_warning("WARNING !!!!!!!!!!!!!!!" // NAMECURRENT // "subroutine obsolete !!!!!!!!!")
!remove!   call log_warning("TODOmd FIND OUT A SOLUTION TO FIX THIS RAHAT") 
!remove!    select case (mld_order)
!remove!     case (mld_linear)
!remove!       if (allocated(w_params)) deallocate (w_params); allocate (w_params(1 + dim_xdesc, 1))
!remove!     case (mld_quadratic)
!remove!       dim_xdesc_quadratic = 1 + dim_xdesc + dim_xdesc**2
!remove!       if (allocated(w_params)) deallocate (w_params); allocate (w_params(dim_xdesc_quadratic, 1))
!remove!    end select
!remove!    _MLD_END_
!remove! end subroutine md_allocate_snap_params

