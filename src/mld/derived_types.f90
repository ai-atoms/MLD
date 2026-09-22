
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

module derived_types

   use module_kind_variables, ONLY: kind_double

   type descriptor_cell

      ! dimension of the descriptor
      integer  :: dim_desc

      ! desc for energy for each atom
      ! energy(1:dim_desc,1:imm)

      ! energy(d,ia) is the d^th dimension of the descriptor on the i^th atom
      real(kind_double), dimension(:, :), allocatable    :: energy
      ! the patch for kernel KNML extension
      real(kind_double), dimension(:, :), allocatable    :: energy_kernel
      ! the patch for kernel 2b extension
      real(kind_double), dimension(:, :), allocatable    :: energy_k2b
      ! the channels for images ...  (ndim*ndim*nchannels)
      real(kind_double), dimension(:, :), allocatable    :: channel


      logical, dimension(:), allocatable :: amask

      ! For composed / hybrid descriptors ...
      integer  :: dim_desc1
      real(kind_double), dimension(:, :), allocatable    :: energy1
      integer  :: dim_desc2
      real(kind_double), dimension(:, :), allocatable    :: energy2

      ! force(1:dim_desc,1:imm,0:imm_neigh,1:3)
      ! force(d,ia,ja,1) the d component of force descriptor for ia^th atom derived
      ! on (ja ,x). If ja=0 then the derivative is with the atom itself.
      real(kind_double), dimension(:, :, :, :), allocatable    :: force
      ! the patch for kernel KNML extension
      real(kind_double), dimension(:, :, :, :), allocatable    :: force_kernel
      ! the patch for kernel 2body extension
      real(kind_double), dimension(:, :, :, :), allocatable    :: force_k2b

      ! For composed / hybrid descriptors...
      real(kind_double), dimension(:, :, :, :), allocatable    :: force1
      real(kind_double), dimension(:, :, :, :), allocatable    :: force2
      ! pack energy - the descriptor for the energy of some configuration (1:dim_desc)
      real(kind_double), dimension(:), allocatable :: pack_energy_linear, &
         pack_energy_quadratic, &
         pack_energy_polyc, &
         pack_energy_kernel, &
         pack_energy_k2b
      ! pack force - the descriptor for force of some atom (1:dim_desc, 3, imm)
      real(kind_double), dimension(:, :, :), allocatable :: pack_force_linear, &
         pack_force_quadratic, &
         pack_force_polyc, &
         pack_force_kernel, &
         pack_force_k2b
      ! pack stress from the derivatives ...(1:dim_desc,1:6)
      real(kind_double), dimension(:, :), allocatable    :: pack_stress_linear, &
         pack_stress_quadratic, &
         pack_stress_polyc, &
         pack_stress_kernel, &
         pack_stress_k2b

      ! statistical distance with respect:
      !   stat_dist_mcd :    imposed Sigma - Hotelling like estimation
      !   stat_dist_maha     Mahalanobis distance
      !   stat _norm         norm of the descritors
      real(kind_double), dimension(:), allocatable :: stat_dist_mcd, stat_dist_maha, stat_norm, stat_norm_mean, stat_energy
      ! n_neigh(ia) number of neighbours of the ia^th atom
      integer, dimension(:), allocatable     :: n_neigh
      ! number of ghost atoms for forces ...
      integer, dimension(:), allocatable     :: n_neigh_ghost
      ! type (species) of neighbours (ia, ith neighbours) and kind (if it is a replica for the atom 1, 2, 3  ...)
      integer, dimension(:, :), allocatable  :: type_neigh, kind_neigh
      logical, dimension(:, :), allocatable  :: incell
      ! the_same thing but gfor ghost atoms ... type (species) of neighbours (ia, ith neighbours) and kind (if it is a replica for the atom 1, 2, 3  ...)
      integer, dimension(:, :), allocatable  :: kind_neigh_ghost, kind_neigh_proc
   end type descriptor_cell

   type system_state
      character(len=80)    :: filename
      character(len=80)    :: hdf5_filename
      character(len=80)    :: hdf5_grpname
      ! a system state
      ! in which class / klm belongs and the i_order of the file in orginal database
      character(len=2)     :: class
      character(len=3)     :: klm
      character(len=6)     :: cnumber
      ! from which line from database the file configuration comes ...
      integer  :: db_line
      ! the number of the configuration in that line
      integer  :: no_file_in_db_line

      integer  :: i_weights
      ! if this config is used for train or test
      ! train=T for training
      ! train=F for testing
      logical  :: train = .false.
      logical  :: selected = .false.
      ! the forces are reliable T/F
      logical  :: has_force = .false.
      ! the energy is reliable  T/F
      logical  :: has_energy = .false.
      ! the stress is reliable T/F
      logical  :: has_stress = .false.
      ! the spins are reliable T/F
      logical  :: has_atomic_spin = .false.
      ! copy of weights from the databse line.
      real(kind_double)    :: w_e, w_f, w_s, w_e_end, w_f_end, w_s_end
      real(kind_double)    :: w_e_model, w_f_model, w_s_model, w_e_end_model, w_f_end_model, w_s_end_model
      ! if the box is small or big
      logical  :: small = .false.

      ! begin small_box all these quantities are active for small boxes

      ! these values will be active in the case small=.true.
      integer  :: nxCell, nyCell, nzCell
      ! number of neighbours for each atom
      integer, dimension(:), allocatable     :: n_neigh, n_neigh_zbl
      ! on which processor the atom is located ...
      integer, dimension(:), allocatable     :: proc_atom
      integer  :: at_start, at_final
      ! type (species) of neighbours (ia, ith neighbours) and kind (if it is a replica for the atom 1, 2, 3  ...)
      integer, dimension(:, :), allocatable  :: type_neigh, kind_neigh, kind_neigh_zbl
      logical, dimension(:, :), allocatable  :: incell
      ! distance between atoms (ia, ith neighbours)
      real(kind=kind(1.d0)), dimension(:, :), allocatable      :: r_ij
      ! distane projected on each direction  u(ia, ith neighbours, 1:3)= r(ia,1:3) - r(ith niegh, 1:3)
      real(kind=kind(1.d0)), dimension(:, :, :), allocatable   :: u_ij, u_per, u_at
      ! end small_box

      character(len=3)     :: state
      ! number of atoms
      integer  :: nat
      ! number of types of atoms
      integer  :: ntypes
      ! volume of the box
      real(kind=kind(1.d0))      :: volume
      ! types of atoms in poscar
      integer, dimension(:), allocatable     :: itype
      ! type of atoms in DB
      integer, dimension(:), allocatable     :: itype_db
      ! type of atoms in fix global
      integer, dimension(:), allocatable     :: itype_to_global
      ! link the type of the atom with the Z in the periodic table
      real(kind_double), dimension(:), allocatable :: fix_type_poscar_to_periodic
      ! direct cell(1:3, component_j), reciprocal bg_cell(1:3,component_j)
      real(kind_double), dimension(3, 3)     :: cell, bg_cell
      real(kind_double), dimension(:, :), allocatable    :: pos_cart, pos_crst
      real(kind_double), dimension(:, :), allocatable    :: force
      real(kind_double), dimension(:, :), allocatable    :: atomic_spin
      ! the reference energy ...
      real(kind_double)    :: ref_energy
      real(kind_double), dimension(:), allocatable :: ref_energy_per_element
      real(kind_double), dimension(3)  :: energy
      real(kind_double), dimension(6)  :: stress
      ! this will be allocated of dimension(ntypes) and contains the atomic mass, Z and covalent radius per type
      real(kind_double), dimension(:), allocatable :: mass_per_type, Z_per_type, covalent_radius_per_type, &
                                                      weight_per_type, weight_per_type_3ch
      logical, dimension(:), allocatable     :: invisible_per_type
      integer  :: spin
      integer  :: im, imm
      integer  :: max_n_neigh = 0

      real(kind_double), dimension(:,:), allocatable :: wspecies_per_type_ch


      real(kind_double), dimension(:), allocatable :: elocal_zbl
      real(kind_double), dimension(:,:,:), allocatable :: flocal_zbl
      real(kind_double), dimension(:,:), allocatable :: slocal_zbl

      real(kind_double) :: ezbl
      real(kind_double), dimension(:,:), allocatable :: fzbl
      real(kind_double), dimension(:), allocatable :: szbl
   end type system_state

   type db_line
      ! an object that englobe one line in the db
      character(len=2)     :: class
      character(len=3)     :: klm
      integer  :: no_total
      integer  :: no_selec
      integer  :: no_start
      ! The T/F on db line that correspond to the energy/force/stress optimization
      logical  :: has_db_force, has_db_energy, has_db_stress
      real(kind_double)    :: w_e, w_f, w_s, w_e_end, w_f_end, w_s_end
      real(kind_double)    :: w_e_model, w_f_model, w_s_model, w_e_end_model, w_f_end_model, w_s_end_model
      ! The T/F on db line that correspond to the energy/force/stress weights optimization
      ! not necessary the same ... e.g. if all configurations on the line have  the flag F there is no optimization for the line
      logical  :: has_optimize_weights_energy, has_optimize_weights_force, has_optimize_weights_stress
      ! maps a compacted 1..no_total index (post drop_short_dist filtering) to
      ! the original on-disk POSCAR file number; unallocated when filtering
      ! is disabled or did not drop anything for this class
      integer, allocatable :: kept_file_number(:)
   end type db_line



   type type_species
      real(kind_double) :: Z1, Z2
      integer :: type1, type2
      character(len=2) :: ch1, ch2
   end type type_species


   type(descriptor_cell), dimension(:), allocatable   :: config_desc, cofig_desc_copy
   type(system_state), dimension(:), allocatable      :: config_real, config_real_copy
   type(db_line), dimension(:), allocatable     :: db_model
   type(type_species), dimension(:), allocatable :: typ_species_half, typ_species_full

end module derived_types
