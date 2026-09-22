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

module module_dump_md_snapshots_extxyz
  use, intrinsic :: iso_fortran_env, dp=>real64
  !use mab_in_ndm_module, only: nlangevin
  integer :: unit_extxyz, unit_energy, unit_msd
  character(len=10) :: atom_type
  real(dp), dimension(:,:), allocatable  :: md_energy
  !real(dp), dimension(:,:), allocatable  :: md_msd
  

  contains

  subroutine init_md_snapshots_extxyz(md_steps)
    implicit none
    integer, intent(in) :: md_steps
    ! Open a file to write. Change 'md_snapshots.extxyz' to your desired file name.
    open(newunit=unit_extxyz, file='md_snapshots.extxyz', status='unknown', action='write')
    open(newunit=unit_energy, file='md_energy.dat', status='unknown', action='write')
    !open(newunit=unit_msd, file='md_msd.dat', status='unknown', action='write')
    if (allocated(md_energy)) deallocate(md_energy)
    allocate(md_energy(3, md_steps))  
      
  end subroutine init_md_snapshots_extxyz

  subroutine close_md_snapshots_extxyz(md_steps)
    implicit none
    integer, intent(in) :: md_steps 
    integer :: ii 
    
    close(unit_extxyz)

    do ii = 1, md_steps
      write(unit_energy, '(i9, 3es25.12)') ii, md_energy(1,ii), md_energy(2,ii), md_energy(3,ii) 
      !write(unit_msd, '(i9, f25.12)') ii, md_msd(ii)
    end do
    close(unit_energy)
    !close(unit_msd)
  end subroutine close_md_snapshots_extxyz


subroutine dump_md_snapshots_extxyz(step_md, energy, positions, simulation_box, N_atoms)

    implicit none
    integer, intent(in) :: step_md, N_atoms
    real(dp), intent(in) :: energy
    real(dp), dimension(3, N_atoms), intent(in) :: positions
    real(dp), dimension(3, 3), intent(in) :: simulation_box
    integer :: i
    !character(len=255) :: line  ! Increased length to ensure enough space

    ! Set the atom type if not already set (example given: 'W')
    atom_type = 'W'

    ! Write the number of atoms
    write(unit_extxyz, *) N_atoms

    ! Construct the comment line including the MD step, energy, and simulation box for this step
    !write(line,*) '(A, I6, A, 9(F12.10,1X), A, F12.10)'
    !write(line, '(A, I6, A, 9(F12.10,1X), A, F12.10)') 'Step=', step_md, ' Lattice="', &
    !     simulation_box(1, 1), simulation_box(1, 2), simulation_box(1, 3), &
    !     simulation_box(2, 1), simulation_box(2, 2), simulation_box(2, 3), &
    !     simulation_box(3, 1), simulation_box(3, 2), simulation_box(3, 3), &
    !     ' Properties=species:S:1:pos:R:3 Energy=', energy

    ! Write the constructed line to the file
    write(unit_extxyz, '(A, I6, A, 9(F20.10,1X), A, A,  F20.10, A)') 'Step=', step_md, ' Lattice="', &
         simulation_box(1, 1), simulation_box(1, 2), simulation_box(1, 3), &
         simulation_box(2, 1), simulation_box(2, 2), simulation_box(2, 3), &
         simulation_box(3, 1), simulation_box(3, 2), simulation_box(3, 3), '"',&
         ' Properties=species:S:1:pos:R:3 Energy=', energy, ' pbc="T T T"'

    ! Write positions of each atom
    do i = 1, N_atoms
        write(unit_extxyz, '(A,3F12.6)') atom_type, positions(1, i), positions(2, i), positions(3, i) 
    end do

end subroutine dump_md_snapshots_extxyz


end module module_dump_md_snapshots_extxyz



module  module_input_md 
   use, intrinsic :: iso_fortran_env, dp=>real64
   real(dp) :: temperature   ! temperature in K 
   real(dp) :: dtmd          ! integration md step in fs = 10^{-15}
   integer  :: eq_steps      ! equilibration steps 
   integer  :: md_steps      ! molecular dynamics steps    
end module  module_input_md 


module module_mld_md
   use, intrinsic :: iso_fortran_env, dp=>real64
   use module_neigh_local, only:  type_neigh_ja
   use mld_logger 
  implicit none
  private
  public :: atom_system, mdsim

  type :: atom
    integer :: ja
    character(len=2) :: chtype     ! Atom type (e.g., 'H', 'O', 'C')
    integer :: itype               ! Integer type for internal use (system box specific)
    integer :: itype_db            ! Integer type from database use
    real(dp), dimension(3) :: pos  ! current position
    real(dp), dimension(3) :: ppos ! previous position
    real(dp), dimension(3) :: vel  ! velocity
    real(dp), dimension(3) :: for  ! force
    type(type_neigh_ja) :: neigh_ja   ! neighbor list
  end type atom

  type :: atom_system
    integer :: md_config
    real(dp), dimension(3,3) :: cell    ! box lattice vectors
    real(dp), dimension(3,3) :: bg_cell ! reciprocal box lattice vectors
    real(dp) :: cutoff                  ! cutoff distance of interactions
    type(atom), dimension(:), allocatable :: atoms
    real(dp), dimension(:,:), allocatable :: spos, svel, sfor ! positions, velocities, forces 3N form. 
    real(dp), dimension(6) :: stress
    real(dp), dimension(:), allocatable :: smass              ! masses of atoms
    real(dp) :: energy_pot, energy_kin, energy_tot
    integer :: natoms
  contains
    procedure :: init_atom_system
    procedure :: update_neighbors_system
    !$! procedure :: export_positions_system
    procedure :: deallocate_system 
  end type atom_system

  type, extends(atom_system) :: mdsim
    real(dp) :: temperature
    real(dp) :: dt
    integer ::  eq_steps, md_steps 
  contains
    procedure :: init_mdsim
    procedure :: propagate_verlet
  end type mdsim

contains

  subroutine init_atom_system(this, iconf)
    use derived_types, only: config_real
    class(atom_system), intent(inout) :: this
    integer, intent(in) :: iconf  ! iconf milady configuration
    integer :: ia

    !number of atoms 
    this%natoms = config_real(iconf)%nat 
    if (allocated(this%atoms)) deallocate(this%atoms) 
    allocate(this%atoms(this%natoms))
    this%md_config = iconf

    if(allocated(this%spos)) deallocate(this%spos) ; allocate(this%spos(3, this%natoms))
    if(allocated(this%svel)) deallocate(this%svel) ; allocate(this%svel(3, this%natoms))
    if(allocated(this%sfor)) deallocate(this%sfor) ; allocate(this%sfor(3, this%natoms))
    if(allocated(this%smass)) deallocate(this%smass) ; allocate(this%smass(this%natoms))

    this%svel(:,:) = 0.0_dp
    this%sfor(:,:) = 0.0_dp

    !positions of atoms 
    do ia = 1, this%natoms
      this%atoms(ia)%ja = ia
      this%atoms(ia)%pos(:) = config_real(iconf)%pos_cart(:, ia)
      this%smass(ia) = config_real(iconf)%mass_per_type(config_real(iconf)%itype(ia))
    end do 
    this%spos(:,:) = config_real(iconf)%pos_cart(:,:) 

    !box lattice vectors
    this%cell(:,:) = config_real(iconf)%cell(:,:)
    this%bg_cell(:,:) = config_real(iconf)%bg_cell(:,:)

  end subroutine init_atom_system

  subroutine update_neighbors_system(this)
    use ml_in_ndm_module, only: imm_neigh
    class(atom_system), intent(inout) :: this
    integer :: i, j, ia_n 
    real(dp) :: dist_sq, cutoff_sq
    real(dp), dimension(3) :: diff 
    real(dp), dimension(1:imm_neigh) :: lr_central
    real(dp), dimension(3,1:imm_neigh) :: ltmp_dxp 
    integer,  dimension(0:imm_neigh) :: li_type, li_central, li_type_db 
    real(dp) :: rr 

    cutoff_sq = this%cutoff**2

    do i = 1, this%natoms
      ia_n = 0
      do j = 1, this%natoms
        if (i /= j) then
          diff = this%atoms(i)%pos - this%atoms(j)%pos
          dist_sq = sum(diff**2)

          if (dist_sq < cutoff_sq) then
            ia_n = ia_n + 1
            ! Add neighbor j to the neighbor list of atom i
            rr = sqrt(dist_sq)
            lr_central(ia_n) = rr 
            ltmp_dxp(:,ia_n) = diff
            li_type(ia_n) = this%atoms(j)%itype
            li_central(ia_n) = j
            li_type_db(ia_n) = this%atoms(j)%itype_db
          end if
        end if
      end do
      this%atoms(i)%neigh_ja%ja_atom = i 
      this%atoms(i)%neigh_ja%max_neigh_local = ia_n
      if (allocated(this%atoms(i)%neigh_ja%r_central)) deallocate(this%atoms(i)%neigh_ja%r_central)
      allocate(this%atoms(i)%neigh_ja%r_central(ia_n))
      this%atoms(i)%neigh_ja%r_central(1:ia_n) = lr_central(1:ia_n)


      if(allocated(this%atoms(i)%neigh_ja%tmp_dxp)) deallocate(this%atoms(i)%neigh_ja%tmp_dxp)
      allocate(this%atoms(i)%neigh_ja%tmp_dxp(3,ia_n))
      this%atoms(i)%neigh_ja%tmp_dxp(:,1:ia_n) = ltmp_dxp(:,1:ia_n)

      li_type(0) = this%atoms(i)%itype
      if (allocated(this%atoms(i)%neigh_ja%i_type)) deallocate(this%atoms(i)%neigh_ja%i_type)
      allocate(this%atoms(i)%neigh_ja%i_type(0:ia_n))
      this%atoms(i)%neigh_ja%i_type(0:ia_n) = li_type(0:ia_n)

      li_central(0) = i
      if (allocated(this%atoms(i)%neigh_ja%i_central)) deallocate(this%atoms(i)%neigh_ja%i_central)
      allocate(this%atoms(i)%neigh_ja%i_central(0:ia_n))
      this%atoms(i)%neigh_ja%i_central(0:ia_n) = li_central(0:ia_n)

      li_type_db(0) = this%atoms(i)%itype_db
      if (allocated(this%atoms(i)%neigh_ja%i_type_db)) deallocate(this%atoms(i)%neigh_ja%i_type_db)
      allocate(this%atoms(i)%neigh_ja%i_type_db(0:ia_n))
      this%atoms(i)%neigh_ja%i_type_db(0:ia_n) = li_type_db(0:ia_n)

    end do
  end subroutine update_neighbors_system

  !$! subroutine export_positions_system(this, filename)
  !$!   class(atom_system), intent(in) :: this
  !$!   character(len=*), intent(in) :: filename
  !$!   integer :: i, unit
  !$!   unit = 10
  !$!   open(unit=unit, file=filename, status='unknown')
  !$!   do i = 1, this%natoms
  !$!     write(unit, *) this%atoms(i)%ja, this%atoms(i)%chtype, this%atoms(i)%pos
  !$!   end do
  !$!   close(unit)
  !$! end subroutine export_positions_system

  subroutine deallocate_system(this)
    class(atom_system), intent(inout) :: this
    integer :: i

    if (allocated(this%atoms)) then
      do i = 1, this%natoms
        if (allocated(this%atoms(i)%neigh_ja%r_central)) deallocate(this%atoms(i)%neigh_ja%r_central)
        if (allocated(this%atoms(i)%neigh_ja%tmp_dxp))   deallocate(this%atoms(i)%neigh_ja%tmp_dxp)
        if (allocated(this%atoms(i)%neigh_ja%i_type))    deallocate(this%atoms(i)%neigh_ja%i_type)
        if (allocated(this%atoms(i)%neigh_ja%i_central)) deallocate(this%atoms(i)%neigh_ja%i_central)
        if (allocated(this%atoms(i)%neigh_ja%i_type_db)) deallocate(this%atoms(i)%neigh_ja%i_type_db)
      end do
      deallocate(this%atoms)
    end if


  end subroutine deallocate_system

  subroutine init_mdsim(this, dt, eq_steps, md_steps, temperature)
    use mld_mpi, only: mld_rank, comm_mld
    class(mdsim), intent(inout) :: this
    real(dp), intent(in) :: dt, temperature 
    integer, intent(in) ::  eq_steps, md_steps 
    integer :: ii 
    real(dp) :: tmp_ekin

    this%temperature = temperature !temperature in K 
    this%dt = dt                   ! time step in fs 
    this%md_steps = md_steps       !number of MD steps
    this%eq_steps = eq_steps       !number of equilibration steps

    ! Initialize velocities with a random Boltzmann distribution
    this%svel(:,:) = 0.0_dp
    if (mld_rank==0) call init_velocities(this, tmp_ekin)
    call comm_mld%bcast(0, tmp_ekin)
    call log_info('ML: MD_mode velocities initialized with kinetic energy: '//vtoa(tmp_ekin))
    do ii = 1, this%natoms
      call comm_mld%bcast(0, this%svel(:, ii))
    end do


  end subroutine init_mdsim

  subroutine propagate_verlet(this)
    use module_units, only:  AMASS
    use mld_mpi, only: mld_rank 
    use module_dump_md_snapshots_extxyz, only: dump_md_snapshots_extxyz, init_md_snapshots_extxyz, close_md_snapshots_extxyz, md_energy
    class(mdsim), intent(inout) :: this
    integer :: i, j
    real(dp) :: tmp_ekin 

    if (mld_rank==0) call init_md_snapshots_extxyz(this%md_steps)
    do i = 1, this%md_steps
      ! Velocity Verlet integration steps
      do j = 1, this%natoms
        !debug write(*,*) 'dt  ', this%dt 
        !debug write(*,*) 'AMASS', AMASS
        !debug write(*,*) 'smass', this%smass(j)
        !debug write(*,*) 'sfor ', this%sfor(:,j)
        !debug write(*,*) 'svel ', this%svel(:,j)
        !debug !stop 
      
        this%svel(:,j) = this%svel(:,j) + 0.5d0 * this%dt * this%sfor(:,j) / (this%smass(j) * AMASS)
        this%spos(:,j) = this%spos(:,j) + this%dt * this%svel(:,j)
      end do

      call calculate_forces(this)
      
      do j = 1, this%natoms
        this%svel(:,j) = this%svel(:,j) + 0.5d0 * this%dt * this%sfor(:,j) / (this%smass(j) * AMASS)
      end do

      tmp_ekin = 0.0_dp 
      do j = 1, this%natoms
        tmp_ekin = tmp_ekin + sum(this%svel(:,j)**2) * (this%smass(j) *   AMASS) / 2.0_dp 
      end do  
      this%energy_kin = tmp_ekin 
      
      this%energy_tot = tmp_ekin  + this%energy_pot 
      if (mld_rank == 0 ) then
        write(*,'("MD step:  ", i7, 3es22.12 )')  i, tmp_ekin, this%energy_pot,  this%energy_tot 
      end if 
      if (mld_rank==0)  md_energy(:,i) = (/ tmp_ekin, this%energy_pot, this%energy_tot /)

      
      if (mld_rank==0) then 
        if (mod(i, 10) == 0) then
          call dump_md_snapshots_extxyz(i, this%energy_tot, this%spos, this%cell, this%natoms)
        end if
      end if 

    end do
    if (mld_rank==0) call close_md_snapshots_extxyz(this%md_steps)

  end subroutine propagate_verlet

  subroutine calculate_forces(this)
    use derived_types, only: config_real
    use snap, only: fp_snap, ene_snap, stress_snap
    class(mdsim), intent(inout) :: this

    config_real(this%md_config)%pos_cart(1:3,1:this%natoms) = this%spos(1:3,1:this%natoms)
    call calfo_ml_as_mld_md
    this%sfor(1:3,1:this%natoms) = fp_snap(1:3,1:this%natoms)
    this%stress(1:6) = stress_snap(1:6)
    this%energy_pot  = ene_snap
  end subroutine calculate_forces

  subroutine init_velocities(this, ekin)
    ! Initialize velocities with random numbers following Boltzmann distribution
    use module_units, only: K_TO_eV, AMASS
    use math, only: sample_1D_gaussian_scalar_box_muller
    class(mdsim), intent(inout) :: this
    real(dp), intent(out)  :: ekin 
    integer :: i, j 
    real(dp) :: kb_temp, stddev, scatmp, tmp, vtmp(3)  
    real(dp) :: mean, low_limit 

    kb_temp = K_TO_eV * this%temperature

    call random_seed()
    do i = 1, this%natoms
      !mean = 0.0_dp ! mean velocity
      stddev = sqrt(kb_temp / ( this%smass(i) * AMASS )) 
      do j =1,3 
      call sample_1D_gaussian_scalar_box_muller(scatmp, mean, stddev)
        vtmp(j) = scatmp 
      end do 
      this%svel(1:3,i) = vtmp(1:3) 
    end do

    low_limit = 5.0_dp * epsilon(1.0_dp)
    do i = 1, this%natoms
      do j = 1, 3
        if (abs(this%svel(j,i)) < low_limit) then
          this%svel(j,i) = 0.0_dp 
        end if
      end do
    end do

    tmp = 0.0_dp
    do i =1, this%natoms
      tmp  = tmp + (this%svel(1,i)**2 + this%svel(2,i)**2 + this%svel(3,i)**2)*this%smass(i) * AMASS/2.0_dp
    end do
    ekin = tmp 

  end subroutine init_velocities

end module module_mld_md

subroutine read_poscar_as_mld_md 
  use, intrinsic :: iso_fortran_env, dp=>real64
  use mld_logger
  use mld_mpi, only: comm_mld, mld_mpi_abort, mld_critical_abort
  use my_mpi_subroutines, only : subworlds_allreduce_int, subworlds_allreduce_vect_int
#ifdef MLD_NDM
  use gen_com_m, only: A2cm, dmtype, fnam
  use gen_com_m_ml, only: im, imm
  use var_pot, only: ntyp
  use ondm_gen_com_m, only: volu
#else
  use ondm_gen_com_m, only: im, imm, fnam
  use ondm_var_pot, only: ntyp
#endif

  use module_db_setup, only: md_iconf, iconf_data, iconf_data_train, iconf_data_test
  use ml_in_ndm_module, only: desc_forces
  use derived_types, only: config_real, config_desc 
  use module_db_poscar, only: fix_ref_energy_per_element
  use module_chemical_species, only: fix_no_of_elements, fix_ch_elements
  use module_md_mld, only: md_allocate_mld_desc
  use snap, only: fp_snap
  use ondm_tab_imm_m, ONLY: dealloc_all_tab_imm,alloc_all_tab_imm

  implicit none 
  integer  :: ii, nb_elements, ij, i_p, n_unmatched
  character(len=:), allocatable    :: fnamposcar 
  logical :: poscarok 


  !initialisation ... very general
  iconf_data = 1
  iconf_data_train = 1
  iconf_data_test = 0

  ! config_real and condif_desc objects
  if (allocated(config_real)) deallocate (config_real); allocate (config_real(iconf_data))
  if (allocated(config_desc)) deallocate (config_desc); allocate (config_desc(iconf_data))
  
  config_real(:)%nat = 0
  md_iconf = 1

  fnamposcar = fnam//'.poscar'
  inquire (file=fnamposcar, exist=poscarok)
  if (poscarok) then 
     config_real(md_iconf)%filename = fnamposcar 
     call log_info('ML:  MD_mode poscar file found: '//fnamposcar)
  else 
     call log_critical('ML:  MD_mode poscar file not found. Please provide a poscar file as input')
     call log_critical('ML:  We expect a file with the name: '//fnamposcar)
     call mld_critical_abort('ML:  MD_mode ...  please provide a poscar file as input')
  end if    

  desc_forces = .false.
  call read_poscar_sasha(md_iconf)
  desc_forces = .true.

  ! Validate: all POSCAR species must exist in the ML model species list.
  ! If a POSCAR species is not in fix_ch_elements, itype_db will be 0 for those
  ! atoms, causing out-of-bounds access in radial spline lookups.
  nb_elements = config_real(md_iconf)%ntypes
  n_unmatched = 0
  do ij = 1, nb_elements
    if (config_real(md_iconf)%itype_to_global(ij) < 1) then
      n_unmatched = n_unmatched + 1
      call log_critical('ML: MD_mode POSCAR species #'//vtoa(ij)// &
                        ' not found in ML model species list (fix_no_of_elements='// &
                        vtoa(fix_no_of_elements)//')')
    end if
  end do
  if (n_unmatched > 0) then
    call log_critical('ML: MD_mode FATAL: '//vtoa(n_unmatched)//' POSCAR species do not match ML model.')
    call log_critical('ML: MD_mode The POSCAR file must contain only species defined in the .ml input file.')
    call log_critical('ML: MD_mode ML model species: '//vtoa(fix_no_of_elements)// &
                      ' elements: '//trim(fix_ch_elements(1)))
    call mld_critical_abort('ML: MD_mode species mismatch between POSCAR and ML model')
  end if

  config_real(md_iconf)%has_energy = .true.
  config_real(md_iconf)%has_force = .true.
  config_real(md_iconf)%has_stress = .true.
  config_real(md_iconf)%force = 0.0_dp 
  config_real(md_iconf)%stress = 0.0_dp

  ! Initialize ref_energy_per_element from fix_ref_energy_per_element
  nb_elements = config_real(md_iconf)%ntypes
  if (allocated(config_real(md_iconf)%ref_energy_per_element)) deallocate(config_real(md_iconf)%ref_energy_per_element)
  allocate(config_real(md_iconf)%ref_energy_per_element(nb_elements))
  config_real(md_iconf)%ref_energy_per_element(:) = 0.0_dp
  if (allocated(fix_ref_energy_per_element)) then
    do ii = 1, min(nb_elements, size(fix_ref_energy_per_element))
      config_real(md_iconf)%ref_energy_per_element(ii) = fix_ref_energy_per_element(ii)
    end do
  end if
  config_real(md_iconf)%ref_energy = 0.0_dp

  call log_info('ML: MD_mode read POSCAR file ... done')

  call fix_poscar_weights
  call fix_atoms_weights(md_iconf)

  ! compatibilities with NDM world ... 
  ntyp = config_real(md_iconf)%ntypes
  call alloc_typ_ml
  call dealloc_all_tab_imm
  call alloc_all_tab_imm(im)
  if (imm .gt. im) then
    call log_warning('imm should be resized to im. MILADY cannot work otherwise '//vtoa([im, imm]))
  end if
  ! compatibilities ...

  call comm_mld%barrier
  ! everyone will need this info during train error
  !$! allocate(tmp_nat(iconf_data))
  !$! tmp_nat(:) = config_real(:)%nat
  !$! call subworlds_allreduce_vect_int(tmp_nat)
  !$! config_real(:)%nat = tmp_nat(:)
  !$! deallocate(tmp_nat)
 
  !debug! write(*,* ) mld_rank, config_real(md_iconf)%nat
  
  call md_allocate_mld_desc(fp_snap)
  !config_real(md_iconf)%volume = volu/A2cm**3
  !debug! write(*,*) mld_rank, config_real(md_iconf)%volume 
  !debug! stop 'oooo'

end subroutine read_poscar_as_mld_md 


subroutine init_subworld_as_mld_md
  use mld_logger
  use mld_mpi, only: mld_rank, mld_size
  use mld_subworld, only: subrank, id_subworld, init_subworld
  use module_db_poscar, only: i_start_conf, i_final_conf, procs_per_file
  use module_db_setup, only: iconf_data
  use set_limits, only: set_limit_for_configs_with_MPI_grid
  implicit none 

  if (procs_per_file /= mld_size ) then
     call log_warning("procs_per_file is different from max rank of MPI world.  For MLD_MD mode this should equal. We force its value to "//vtoa(mld_size)) 
     procs_per_file = mld_size
  end if 
  call init_subworld(procs_per_file)
  call set_limit_for_configs_with_MPI_grid(iconf_data, i_start_conf, i_final_conf)
  write(*,'("ML grid info w-rank, sub-rank, id-sub, i_start_conf, ifinal_conf :", 5i6)') mld_rank, subrank, id_subworld, i_start_conf, i_final_conf
  !train_time=.false.
end subroutine init_subworld_as_mld_md

subroutine calfo_ml_as_mld_md
  use, intrinsic :: iso_fortran_env, dp=>real64
  use mld_logger 
  use module_db_poscar, only: i_start_conf, i_final_conf
  use mld_mpi, only: comm_mld
  use time_check_general, only: test_time_total, test_time_desc, &
         test_time_eval, test_time_neigh, test_tot_time, MY_MPI_WTIME, &
         debug_time, tot_time
  use main_mld_mod, only: md_mld_compute_energy, md_mld_compute_force, md_mld_compute_stress 
  use derived_types, only: config_real    
  use module_md_mld, only: md_allocate_mld_desc 
  use module_kernel_2b, only: test_time_for_desc_kernel_2b, time_for_desc_kernel_2b, activate_k2b
  use ml_in_ndm_module, only:  ml_type, ml_type_krr
  use module_kernel, only: time_for_desc_kernel, test_time_for_desc_kernel
  use snap, only: fp_snap
  implicit none
  real(dp) :: t0, t1, t2, t3, t4, t5
  integer :: i
  real(dp) :: local_ref_energy, local_ref_zbl
  logical :: post_desc


  _NAMECURRENT_("calfo_ml_as_mld_md")
  _MLD_BEGIN_

  t0 = MY_MPI_WTIME()

  !$! if (force_variance) call compute_variance_matrix
      
  do i = i_start_conf, i_final_conf
    t1 = MY_MPI_WTIME()
    call test_if_config_is_small(i)
    call calc_neighbours(i)
    t2 = MY_MPI_WTIME()
    test_time_neigh = test_time_neigh + (t2 - t1)
    !!!write (6,*) 'before md_allocate', i, rangml
    call md_allocate_mld_desc(fp_snap)
    !!!write (6,*) 'after md_allocate', i, rangml
    post_desc=.true.
    call compute_descriptors(i, post_desc)
    t3 = MY_MPI_WTIME()
    test_time_desc = test_time_desc + (t3 - t2)
    !debug if (rangml==0) write (6,*) 'after',
    !mpi_rangml =0
    !mpi_rangml
    !if (rangml==0) then
    call md_mld_compute_energy(i, pack_opt=.true.)
    local_ref_energy = config_real(i)%ref_energy
    local_ref_zbl = config_real(i)%ezbl
    !y_e_test_base(i_e_test_snap) = config_real(i)%energy(iread_energy) + local_ref_energy
    call md_mld_compute_force(i, pack_opt=.true.)
    call md_mld_compute_stress(i, pack_opt=.true.)
    call train_deallocate_desc(i)
    t4 = MY_MPI_WTIME()
    test_time_eval = test_time_eval + (t4 - t3)
    end do

    call comm_mld%barrier()
    t5 = MY_MPI_WTIME()
    test_time_total = test_time_total + (t5 - t0)
    if (ml_type == ml_type_krr) then
       test_time_for_desc_kernel = time_for_desc_kernel
    end if
    if (activate_k2b) then
       test_time_for_desc_kernel_2b = time_for_desc_kernel_2b
    end if

    if (debug_time) then
       test_tot_time = tot_time
    end if

    _MLD_END_
end subroutine calfo_ml_as_mld_md

subroutine repport_time_mld_md 
  use ml_in_ndm_module, only: rangml, ml_type, ml_type_krr,  train_time
  use time_check_general, only: test_time_total, &
                                test_time_desc, &
                                test_time_eval, &
                                test_time_neigh, debug_time, test_tot_time
  use module_kernel, only: test_time_for_desc_kernel
  use module_kernel_2b, only: activate_k2b, test_time_for_desc_kernel_2b                    
  implicit none

  if (rangml == 0) write (6, '("ML:------------------REPORTING TIME---------------------------")')
  
  call repport_time(2, 0.d0, test_time_total, "ML: TEST TIME")
  call repport_time(3, 0.d0, test_time_neigh, "ML: time for searching neighbours ")
  call repport_time(3, 0.d0, test_time_desc, "ML: time for computing desc")
  if (debug_time) then 
    train_time = .false.   
    call repport_debug_time_descriptor(test_tot_time)
  end if 
      
  if (ml_type == ml_type_krr) then
    call repport_time(4, 0.d0, test_time_desc - test_time_for_desc_kernel, "ML: pure desc")
    call repport_time(4, 0.d0, test_time_for_desc_kernel, "ML: kernel desc")
  end if
  if (activate_k2b)  call repport_time(3, 0.d0, test_time_for_desc_kernel_2b, "ML: kernel 2-body")
  if (debug_time) call repport_debug_time_kernel()
  
  call repport_time(3, 0.d0, test_time_eval, "ML: time for evaluation ")
end subroutine repport_time_mld_md


subroutine rrrroule_ma_poule
!$! program molecular_dynamics
  use mld_logger
  use module_mld_md, only: mdsim
  use module_db_setup, only: md_iconf
  use mld_subworld, only: close_subworld
  use module_input_md, only: md_steps, eq_steps, dtmd, temperature

  implicit none
  type(mdsim) :: mdrun

  call read_poscar_as_mld_md
  call init_subworld_as_mld_md
  !call calfo_ml_as_mld_md
  call mdrun%init_atom_system(md_iconf)
  call mdrun%init_mdsim(dtmd, eq_steps, md_steps, temperature)
  call mdrun%propagate_verlet
  call repport_time_mld_md
  call close_subworld()

!$! 
!$!   ! Deallocate system
!$!   call integrator%deallocate_system()
end subroutine rrrroule_ma_poule


subroutine read_mld_md 
  use mld_logger
  use module_input_md, only: dtmd, temperature, eq_steps, md_steps
  use ml_in_ndm_module, only: mld_dmtype, MD_MLD_DMTYPE, ML_MLD_DMTYPE ! 181 and 18 
#ifdef MLD_NDM
   use gen_com_m, ONLY: fnam,  lenfnam
#else
   use ondm_gen_com_m, ONLY: fnam
#endif
  implicit none 
  character(len=:), allocatable    :: fnamtin
  integer :: lumdml 
  logical :: mdmlok

  !TODOmd MD_mode read in rnak 0 and then broadcast ... 
  namelist / input_mdml / eq_steps, md_steps, temperature, dtmd
#ifdef MLD_NDM
  fnamtin = fnam(1:lenfnam)//'.mdml'
#else
  fnamtin = fnam//'.mdml'
#endif

  mdmlok = .false.
  inquire (file=fnamtin, exist=mdmlok)
  if (mdmlok) then 
    eq_steps = 100
    md_steps = 1000
    temperature = 300.0d0
    dtmd = 0.1d0

    mld_dmtype = MD_MLD_DMTYPE
    open (newunit=lumdml, file=fnamtin, status='unknown')
    read (lumdml, nml=input_mdml)
    close (lumdml)

    call log_info('ML: MDML input file ... ' // trim(fnamtin) // ' ... read')

  else
    mld_dmtype = ML_MLD_DMTYPE 
    return 
  end if 

  call log_info("ML: read input internal mld_dmtype ....:"//vtoa(mld_dmtype))

  if (mld_dmtype == MD_MLD_DMTYPE) then 
     call log_info(" ------------------------------------------------ ")
     call log_info("|                                                |")
     call log_info("|>>>>>         MILADY in MD MODE            <<<<<|")
     call log_info("|                                                |")
     call log_info(" ------------------------------------------------ ")
  end if 

end subroutine read_mld_md 