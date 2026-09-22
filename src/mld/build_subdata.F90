! HND XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
! HND X
! HND X   MiLaDy (Machine Learning Dynamics)
! HND X
! HND X
! HND X   Milady was written and designed by Mihai-Cosmin Marinica, Alexandra M. Goryaeva
! HND X   Copyright 2015-2024.
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

! ============================================================================
! MODULE module_extxyz
!   Extended XYZ database reader.
!   Allows reading databases directly from .xyz / .extxyz files
!   instead of the traditional POSCAR folder + db_model.in format.
!   Activated when db_path ends with .xyz or .extxyz
! ============================================================================
module module_extxyz
   use module_kind_variables, only: kind_double
   implicit none
   private

   !> Whether the database is in extended XYZ format (auto-detected from db_path)
   logical, public :: db_xyz = .false.

   !> Maximum length of info line
   integer, parameter :: MAX_INFO_LEN = 4096

   !> Per-config metadata stored after the pre-scan
   type, public :: xyz_config_info
      integer            :: natoms       ! number of atoms in this config
      character(len=2)   :: class_id     ! class label (from "class=" key, default "01")
      character(len=6)   :: cnumber      ! 6-digit config number within class
      logical            :: has_energy   ! energy= key found
      logical            :: has_forces   ! forces in Properties
      logical            :: has_stress   ! stress= key found
      logical            :: has_lattice  ! Lattice= key found
      real(kind_double)  :: energy       ! total energy
      real(kind_double)  :: free_energy  ! free energy (if present, else = energy)
      real(kind_double)  :: stress(9)    ! 3x3 stress in row-major (xx,xy,xz,yx,yy,yz,zx,zy,zz)
      real(kind_double)  :: lattice(3,3) ! cell vectors as columns: lattice(:,1)=a, lattice(:,2)=b, lattice(:,3)=c
      integer            :: file_line    ! line number of the natoms line (1-based) for this config
      integer(8)         :: byte_offset  ! byte position of the first atom line (for stream I/O)
      logical            :: keep = .true. ! .false. if dropped by drop_short_dist filtering
   end type xyz_config_info

   !> Array of pre-scanned config metadata (allocated in xyz_scan_database)
   type(xyz_config_info), dimension(:), allocatable, public :: xyz_configs

   !> Total number of configs found in the XYZ file
   integer, public :: xyz_n_configs_total = 0

   !> Mapping: xyz_config_map(i) = index into xyz_configs for the i-th config_real slot
   integer, dimension(:), allocatable, public :: xyz_config_map

   !> First-appearance-order auto-registration of subset= names (used when a
   !> config has no class= key). Dataset-independent: no fixed name table,
   !> class ids '01'..'99' are assigned to whatever names are actually found.
   !> Reset at the start of every xyz_scan_database() call.
   integer :: n_auto_subset = 0
   character(len=64) :: auto_subset_name(99)
   character(len=2)  :: auto_subset_class(99)

   public :: xyz_detect_mode
   public :: xyz_scan_database
   public :: xyz_prepare_database
   public :: read_xyz_config

contains

   !---------------------------------------------------------------------------
   !> Detect whether db_path points to an XYZ file
   !---------------------------------------------------------------------------
   subroutine xyz_detect_mode()
      use module_db_setup, only: db_path
      implicit none
      integer :: l

      db_xyz = .false.
      if (.not. allocated(db_path)) return
      l = len_trim(db_path)
      if (l < 4) return

      ! Check for .extxyz extension (7 chars)
      if (l >= 7) then
         if (db_path(l-6:l) == '.extxyz') then
            db_xyz = .true.
            return
         end if
      end if
      ! Check for .xyz extension (4 chars)
      if (db_path(l-3:l) == '.xyz') then
         db_xyz = .true.
         return
      end if
   end subroutine xyz_detect_mode


   !---------------------------------------------------------------------------
   !> Pre-scan the XYZ file to count configs, extract metadata from info lines.
   !> Only rangml==0 reads; results are broadcast via MPI_BCAST on mpi_comm_mld.
   !---------------------------------------------------------------------------
   subroutine xyz_scan_database()
      use module_db_setup, only: db_path, drop_short_dist
      use mld_mpi, only: mpi_comm_mld, rangml, mld_ierror
      use mpi
      use mld_logger
      implicit none

      integer :: inp, ios, natoms, iconf, i, line_number
      character(len=MAX_INFO_LEN) :: info_line
      character(len=500) :: atom_line
      character(len=20) :: ctmp
      character(len=2) :: sp_dummy
      integer :: n_total
      logical :: lexist
      ! Class counting
      integer :: class_counts(100)
      character(len=2) :: class_labels(100)
      integer :: n_classes, ic
      ! drop_short_dist filtering
      real(kind_double), allocatable :: scan_pos(:,:)
      real(kind_double) :: scan_box(3,3)
      logical :: too_close
      integer :: n_dropped

      _NAMECURRENT_("xyz_scan_database")
      _MLD_BEGIN_

      ! ------- First pass: count configs (rank 0 only) -------
      n_total = 0
      if (rangml == 0) then
         inquire(file=trim(db_path), exist=lexist)
         if (.not. lexist) then
            write(6,*) 'ML: FATAL: XYZ file not found: ', trim(db_path)
            call MPI_Abort(mpi_comm_mld, 1, mld_ierror)
         end if
         open(newunit=inp, file=trim(db_path), status='old', action='read', iostat=ios)
         if (ios /= 0) then
            write(6,*) 'ML: FATAL: cannot open XYZ file: ', trim(db_path)
            call MPI_Abort(mpi_comm_mld, 1, mld_ierror)
         end if

         line_number = 0
         do
            read(inp, '(a)', iostat=ios) ctmp
            if (ios /= 0) exit
            line_number = line_number + 1
            ctmp = adjustl(ctmp)
            read(ctmp, *, iostat=ios) natoms
            if (ios /= 0) cycle
            if (natoms <= 0) cycle
            ! Skip info line
            read(inp, '(a)', iostat=ios) info_line
            if (ios /= 0) exit
            line_number = line_number + 1
            ! Skip atom lines
            do i = 1, natoms
               read(inp, '(a)', iostat=ios) info_line
               if (ios /= 0) exit
               line_number = line_number + 1
            end do
            if (ios /= 0) exit
            n_total = n_total + 1
         end do
         close(inp)
         write(6, '("ML: XYZ database scan: found ", i0, " configurations in ", a)') &
            n_total, trim(db_path)
      end if

      ! Broadcast n_total to all ranks
      call MPI_BCAST(n_total, 1, MPI_INTEGER, 0, mpi_comm_mld, mld_ierror)
      xyz_n_configs_total = n_total

      if (n_total == 0) then
         if (rangml == 0) write(6,*) 'ML: FATAL: no configurations found in XYZ file'
         call MPI_Abort(mpi_comm_mld, 1, mld_ierror)
      end if

      ! Allocate config info array on all ranks
      if (allocated(xyz_configs)) deallocate(xyz_configs)
      allocate(xyz_configs(n_total))

      ! ------- Second pass: extract metadata (rank 0 only) -------
      if (rangml == 0) then
         open(newunit=inp, file=trim(db_path), status='old', action='read', access='stream', form='formatted', iostat=ios)
         n_classes = 0
         class_counts(:) = 0
         class_labels(:) = '  '
         iconf = 0
         line_number = 0
         n_dropped = 0
         n_auto_subset = 0
         auto_subset_name(:) = ' '
         auto_subset_class(:) = '  '

         do
            read(inp, '(a)', iostat=ios) ctmp
            if (ios /= 0) exit
            line_number = line_number + 1
            ctmp = adjustl(ctmp)
            read(ctmp, *, iostat=ios) natoms
            if (ios /= 0) cycle
            if (natoms <= 0) cycle

            iconf = iconf + 1
            xyz_configs(iconf)%natoms = natoms
            xyz_configs(iconf)%file_line = line_number
            xyz_configs(iconf)%keep = .true.

            ! Read info line
            read(inp, '(a)', iostat=ios) info_line
            if (ios /= 0) exit
            line_number = line_number + 1

            ! Parse the info line
            call parse_xyz_info_line(info_line, xyz_configs(iconf))

            ! Record byte position of first atom line (portable Fortran 2003 stream)
            inquire(unit=inp, pos=xyz_configs(iconf)%byte_offset)

            ! Register the class label now (even if this instance ends up
            ! dropped below), so a class only ever populated by dropped
            ! configs still shows up in the summary with a zero count.
            ic = find_or_add_class(xyz_configs(iconf)%class_id, &
                                   class_labels, class_counts, n_classes)

            ! Atom lines: parsed for the drop_short_dist check when active,
            ! otherwise just skipped (cheap default path, unchanged behavior).
            if (drop_short_dist > 0.d0) then
               if (allocated(scan_pos)) deallocate(scan_pos)
               allocate(scan_pos(3, natoms))
               do i = 1, natoms
                  read(inp, '(a)', iostat=ios) atom_line
                  if (ios /= 0) exit
                  line_number = line_number + 1
                  read(atom_line, *, iostat=ios) sp_dummy, scan_pos(1,i), scan_pos(2,i), scan_pos(3,i)
                  if (ios /= 0) scan_pos(:,i) = 0.d0
               end do
               if (ios == 0) then
                  ! has_lattice alone is the correct test (see read_xyz_config
                  ! for the same fix): a diagonal-magnitude check misfires on
                  ! valid zero-diagonal cells such as the FCC primitive-cell
                  ! convention a=(0,L,L), b=(L,0,L), c=(L,L,0).
                  if (xyz_configs(iconf)%has_lattice) then
                     scan_box = xyz_configs(iconf)%lattice
                  else
                     scan_box = 0.d0
                     scan_box(1,1) = 100.d0; scan_box(2,2) = 100.d0; scan_box(3,3) = 100.d0
                  end if
                  call config_has_short_pair(natoms, scan_pos, scan_box, drop_short_dist, too_close)
                  if (too_close) then
                     xyz_configs(iconf)%keep = .false.
                     n_dropped = n_dropped + 1
                  end if
               end if
            else
               ! Skip atom lines (no filtering requested)
               do i = 1, natoms
                  read(inp, '(a)', iostat=ios) info_line
                  if (ios /= 0) exit
                  line_number = line_number + 1
               end do
            end if
            if (ios /= 0) exit

            ! Only count/number kept configs, so this stays contiguous
            ! (no gaps), matching the class_list built later in
            ! xyz_prepare_database.
            if (xyz_configs(iconf)%keep) then
               class_counts(ic) = class_counts(ic) + 1
               write(xyz_configs(iconf)%cnumber, '(i6.6)') class_counts(ic)
            end if
         end do
         close(inp)
         if (allocated(scan_pos)) deallocate(scan_pos)

         write(6, '("ML: XYZ scan complete: ", i0, " configs, ", i0, " classes")') &
            iconf, n_classes
         do ic = 1, n_classes
            write(6, '("ML:   class ", a2, ": ", i0, " configs")') &
               class_labels(ic), class_counts(ic)
         end do
         if (drop_short_dist > 0.d0) then
            write(6, '("ML: drop_short_dist = ", f10.5, " Angstrom: dropped ", i0, &
               &" of ", i0, " XYZ configs (interatomic distance below threshold)")') &
               drop_short_dist, n_dropped, iconf
         end if
      end if

      ! Broadcast all config metadata to all ranks
      if (rangml == 0) write(6,'("ML: XYZ DEBUG: before bcast_all_configs")')
      call flush(6)
      call xyz_bcast_all_configs(n_total)
      if (rangml == 0) write(6,'("ML: XYZ DEBUG: after bcast_all_configs")')
      call flush(6)

      _MLD_END_
   end subroutine xyz_scan_database


   !---------------------------------------------------------------------------
   !> Broadcast xyz_configs array from rank 0 to all ranks using MPI_BCAST
   !---------------------------------------------------------------------------
   subroutine xyz_bcast_all_configs(n)
      use mld_mpi, only: mpi_comm_mld, mld_ierror
      use mpi
      implicit none
      integer, intent(in) :: n
      integer :: i, idx
      integer, allocatable :: ibuf(:), lbuf_int(:)
      real(kind_double), allocatable :: rbuf(:)
      ! Use integer arrays to broadcast class_id (2 chars) and cnumber (6 chars)
      ! Avoids issues with deferred-length allocatable character + MPI_BCAST
      integer, allocatable :: cbuf_class(:), cbuf_cnum(:)

      ! --- Integers: natoms and file_line ---
      allocate(ibuf(2*n))
      do i = 1, n
         ibuf(i)   = xyz_configs(i)%natoms
         ibuf(n+i) = xyz_configs(i)%file_line
      end do
      call MPI_BCAST(ibuf, 2*n, MPI_INTEGER, 0, mpi_comm_mld, mld_ierror)
      do i = 1, n
         xyz_configs(i)%natoms    = ibuf(i)
         xyz_configs(i)%file_line = ibuf(n+i)
      end do
      deallocate(ibuf)

      ! --- Byte offsets (integer(8), broadcast as MPI_INTEGER8) ---
      block
         integer(8), allocatable :: i8buf(:)
         allocate(i8buf(n))
         do i = 1, n
            i8buf(i) = xyz_configs(i)%byte_offset
         end do
         call MPI_BCAST(i8buf, n, MPI_INTEGER8, 0, mpi_comm_mld, mld_ierror)
         do i = 1, n
            xyz_configs(i)%byte_offset = i8buf(i)
         end do
         deallocate(i8buf)
      end block

      ! --- Logicals (as integers: 1=true, 0=false) ---
      allocate(lbuf_int(5*n))
      do i = 1, n
         lbuf_int(i)       = merge(1, 0, xyz_configs(i)%has_energy)
         lbuf_int(n+i)     = merge(1, 0, xyz_configs(i)%has_forces)
         lbuf_int(2*n+i)   = merge(1, 0, xyz_configs(i)%has_stress)
         lbuf_int(3*n+i)   = merge(1, 0, xyz_configs(i)%has_lattice)
         lbuf_int(4*n+i)   = merge(1, 0, xyz_configs(i)%keep)
      end do
      call MPI_BCAST(lbuf_int, 5*n, MPI_INTEGER, 0, mpi_comm_mld, mld_ierror)
      do i = 1, n
         xyz_configs(i)%has_energy  = (lbuf_int(i)       == 1)
         xyz_configs(i)%has_forces  = (lbuf_int(n+i)     == 1)
         xyz_configs(i)%has_stress  = (lbuf_int(2*n+i)   == 1)
         xyz_configs(i)%has_lattice = (lbuf_int(3*n+i)   == 1)
         xyz_configs(i)%keep        = (lbuf_int(4*n+i)   == 1)
      end do
      deallocate(lbuf_int)

      ! --- Reals: energy, free_energy, stress(9), lattice(9) = 20 per config ---
      allocate(rbuf(20*n))
      do i = 1, n
         idx = (i-1)*20
         rbuf(idx+1)  = xyz_configs(i)%energy
         rbuf(idx+2)  = xyz_configs(i)%free_energy
         rbuf(idx+3:idx+11) = xyz_configs(i)%stress(1:9)
         rbuf(idx+12) = xyz_configs(i)%lattice(1,1)
         rbuf(idx+13) = xyz_configs(i)%lattice(2,1)
         rbuf(idx+14) = xyz_configs(i)%lattice(3,1)
         rbuf(idx+15) = xyz_configs(i)%lattice(1,2)
         rbuf(idx+16) = xyz_configs(i)%lattice(2,2)
         rbuf(idx+17) = xyz_configs(i)%lattice(3,2)
         rbuf(idx+18) = xyz_configs(i)%lattice(1,3)
         rbuf(idx+19) = xyz_configs(i)%lattice(2,3)
         rbuf(idx+20) = xyz_configs(i)%lattice(3,3)
      end do
      call MPI_BCAST(rbuf, 20*n, MPI_DOUBLE_PRECISION, 0, mpi_comm_mld, mld_ierror)
      do i = 1, n
         idx = (i-1)*20
         xyz_configs(i)%energy      = rbuf(idx+1)
         xyz_configs(i)%free_energy = rbuf(idx+2)
         xyz_configs(i)%stress(1:9) = rbuf(idx+3:idx+11)
         xyz_configs(i)%lattice(1,1) = rbuf(idx+12)
         xyz_configs(i)%lattice(2,1) = rbuf(idx+13)
         xyz_configs(i)%lattice(3,1) = rbuf(idx+14)
         xyz_configs(i)%lattice(1,2) = rbuf(idx+15)
         xyz_configs(i)%lattice(2,2) = rbuf(idx+16)
         xyz_configs(i)%lattice(3,2) = rbuf(idx+17)
         xyz_configs(i)%lattice(1,3) = rbuf(idx+18)
         xyz_configs(i)%lattice(2,3) = rbuf(idx+19)
         xyz_configs(i)%lattice(3,3) = rbuf(idx+20)
      end do
      deallocate(rbuf)

      ! --- Characters: class_id (2 chars) and cnumber (6 chars) as ichar arrays ---
      allocate(cbuf_class(2*n), cbuf_cnum(6*n))
      cbuf_class = 0
      cbuf_cnum  = 0
      do i = 1, n
         cbuf_class((i-1)*2+1) = ichar(xyz_configs(i)%class_id(1:1))
         cbuf_class((i-1)*2+2) = ichar(xyz_configs(i)%class_id(2:2))
         cbuf_cnum((i-1)*6+1) = ichar(xyz_configs(i)%cnumber(1:1))
         cbuf_cnum((i-1)*6+2) = ichar(xyz_configs(i)%cnumber(2:2))
         cbuf_cnum((i-1)*6+3) = ichar(xyz_configs(i)%cnumber(3:3))
         cbuf_cnum((i-1)*6+4) = ichar(xyz_configs(i)%cnumber(4:4))
         cbuf_cnum((i-1)*6+5) = ichar(xyz_configs(i)%cnumber(5:5))
         cbuf_cnum((i-1)*6+6) = ichar(xyz_configs(i)%cnumber(6:6))
      end do
      call MPI_BCAST(cbuf_class, 2*n, MPI_INTEGER, 0, mpi_comm_mld, mld_ierror)
      call MPI_BCAST(cbuf_cnum,  6*n, MPI_INTEGER, 0, mpi_comm_mld, mld_ierror)
      do i = 1, n
         xyz_configs(i)%class_id(1:1) = char(cbuf_class((i-1)*2+1))
         xyz_configs(i)%class_id(2:2) = char(cbuf_class((i-1)*2+2))
         xyz_configs(i)%cnumber(1:1)  = char(cbuf_cnum((i-1)*6+1))
         xyz_configs(i)%cnumber(2:2)  = char(cbuf_cnum((i-1)*6+2))
         xyz_configs(i)%cnumber(3:3)  = char(cbuf_cnum((i-1)*6+3))
         xyz_configs(i)%cnumber(4:4)  = char(cbuf_cnum((i-1)*6+4))
         xyz_configs(i)%cnumber(5:5)  = char(cbuf_cnum((i-1)*6+5))
         xyz_configs(i)%cnumber(6:6)  = char(cbuf_cnum((i-1)*6+6))
      end do
      deallocate(cbuf_class, cbuf_cnum)

   end subroutine xyz_bcast_all_configs


   !---------------------------------------------------------------------------
   !> Set up the database structures (db_model, config_real, etc.) from XYZ scan.
   !> This replaces prepare_database + prepare_name_file_from_db for XYZ mode.
   !---------------------------------------------------------------------------
   subroutine xyz_prepare_database()
      !! Read db_model.in for train/test selection (same as POSCAR path),
      !! then build a mapping from config_real index → xyz_configs index.
      !! The existing read_db_file / prepare_name_file_from_db are reused
      !! so that selection_type, no_total, no_selec, weights, T/F flags
      !! all work identically to the POSCAR workflow.
      use ml_in_ndm_module, only: rangml, debug
      use module_db_setup, only: db_file, db_path, db_train, db_test, &
         iconf_data, iconf_data_train, iconf_data_test
      use derived_types, only: db_model, config_real, config_desc
      use mld_unit
      use mld_logger
      use mld_mpi
      implicit none

      integer :: i, ic, n_total, n_classes, iunit, ierr, n_db_class
      logical :: lexist
      ! Per-class ordered lists of xyz_configs indices
      integer :: class_count_xyz(100)        ! count of configs per class in XYZ file
      character(len=2) :: class_labels(100)  ! class labels
      integer, allocatable :: class_list(:,:) ! class_list(ic, j) = xyz_configs index of j-th config in class ic
      integer :: max_per_class

      _NAMECURRENT_("xyz_prepare_database")
      _MLD_BEGIN_

      n_total = xyz_n_configs_total

      ! ---- Build per-class ordered lists from xyz_configs ----
      n_classes = 0
      class_count_xyz(:) = 0
      class_labels(:) = '  '
      do i = 1, n_total
         if (.not. xyz_configs(i)%keep) cycle
         ic = 0
         do ic = 1, n_classes
            if (xyz_configs(i)%class_id == class_labels(ic)) exit
         end do
         if (ic > n_classes) then
            n_classes = n_classes + 1
            class_labels(n_classes) = xyz_configs(i)%class_id
            ic = n_classes
         end if
         class_count_xyz(ic) = class_count_xyz(ic) + 1
      end do

      ! Allocate per-class index lists (dropped configs are excluded, so
      ! numbering here stays contiguous/no-gaps for db_model.in selection).
      max_per_class = maxval(class_count_xyz(1:n_classes))
      allocate(class_list(n_classes, max_per_class))
      class_list = 0
      class_count_xyz(:) = 0  ! reset for second pass
      do i = 1, n_total
         if (.not. xyz_configs(i)%keep) cycle
         do ic = 1, n_classes
            if (xyz_configs(i)%class_id == class_labels(ic)) exit
         end do
         class_count_xyz(ic) = class_count_xyz(ic) + 1
         class_list(ic, class_count_xyz(ic)) = i
      end do

      if (rangml == 0) then
         do ic = 1, n_classes
            write(6, '("ML: XYZ class ", a, ": ", i0, " configs in XYZ file")') &
               class_labels(ic), class_count_xyz(ic)
         end do
      end if

      ! ---- Read db_model.in via existing POSCAR routines ----
      inquire(file=db_file, exist=lexist)
      if (.not. lexist) then
         call log_critical('xyz_prepare_database: db_file not found: '//trim(db_file))
         call mld_mpi_finalize('xyz_prepare_database: provide db_model.in')
         stop
      end if
      open(newunit=iunit, file=db_file, status='old', action='read', iostat=ierr)
      if (ierr /= 0) then
         call log_critical('xyz_prepare_database: error opening db_file: '//trim(db_file))
         call mld_mpi_finalize('xyz_prepare_database: check file')
         stop
      end if

      ! This reads db_model.in, allocates db_model, config_real, config_desc,
      ! db_train, db_test, sets iconf_data, iconf_data_train, iconf_data_test
      call read_n_db_class(iunit, n_db_class)
      call read_db_file(iunit, n_db_class)
      close(iunit)

      ! ---- Validate: check that db_model classes exist in XYZ and no_total fits ----
      do i = 1, size(db_model)
         ic = 0
         do ic = 1, n_classes
            if (db_model(i)%class == class_labels(ic)) exit
         end do
         if (ic > n_classes) then
            if (rangml == 0) write(6, '("ML: ERROR: class ", a, " in db_model.in not found in XYZ file")') &
               db_model(i)%class
            call mld_mpi_finalize('xyz_prepare_database: class not found in XYZ')
            stop
         end if
         if (db_model(i)%no_total > class_count_xyz(ic)) then
            if (rangml == 0) then
               write(6, '("ML: WARNING: class ", a, " db_model no_total=", i0, &
                  &" exceeds XYZ count=", i0, ", clamping")') &
                  db_model(i)%class, db_model(i)%no_total, class_count_xyz(ic)
            end if
            db_model(i)%no_total = class_count_xyz(ic)
            if (db_model(i)%no_selec > db_model(i)%no_total) then
               db_model(i)%no_selec = db_model(i)%no_total
            end if
            ! Recompute iconf_data etc. after clamping
         end if
      end do

      ! Recompute iconf_data/train/test after possible clamping
      iconf_data_train = 0
      iconf_data = 0
      do i = 1, size(db_model)
         iconf_data_train = iconf_data_train + db_model(i)%no_selec
         iconf_data = iconf_data + db_model(i)%no_total - db_model(i)%no_start + 1
      end do
      iconf_data_test = iconf_data - iconf_data_train

      ! Re-allocate arrays with corrected sizes
      if (allocated(config_real)) deallocate(config_real); allocate(config_real(iconf_data))
      if (allocated(config_desc)) deallocate(config_desc); allocate(config_desc(iconf_data))
      if (allocated(db_train)) deallocate(db_train); allocate(db_train(iconf_data_train))
      if (allocated(db_test))  deallocate(db_test);  allocate(db_test(max(iconf_data_test,0)))

      ! ---- Use prepare_name_file_from_db to do selection (train/test split) ----
      ! This fills config_real(:)%class, klm, cnumber, filename, train, selected,
      ! w_e, w_f, w_s, has_energy/force/stress, db_train(:), db_test(:)
      call prepare_name_file_from_db()

      ! ---- Build xyz_config_map: config_real(i) → xyz_configs index ----
      ! config_real(i)%class is the class label
      ! config_real(i)%no_file_in_db_line is the 1-based sequence number within that class
      if (allocated(xyz_config_map)) deallocate(xyz_config_map)
      allocate(xyz_config_map(iconf_data))
      do i = 1, iconf_data
         ! Find which class
         do ic = 1, n_classes
            if (config_real(i)%class == class_labels(ic)) exit
         end do
         ! Map sequence number to xyz_configs index
         xyz_config_map(i) = class_list(ic, config_real(i)%no_file_in_db_line)
      end do

      ! Override has_energy/force/stress from actual XYZ data
      do i = 1, iconf_data
         config_real(i)%has_energy = config_real(i)%has_energy .and. &
            xyz_configs(xyz_config_map(i))%has_energy
         config_real(i)%has_force  = config_real(i)%has_force  .and. &
            xyz_configs(xyz_config_map(i))%has_forces
         config_real(i)%has_stress = config_real(i)%has_stress .and. &
            xyz_configs(xyz_config_map(i))%has_stress
      end do

      call get_distinct_classes()

      ! ---- Compact xyz_configs: keep only referenced entries, release the rest ----
      block
         logical, allocatable :: used(:)
         type(xyz_config_info), allocatable :: xyz_configs_compact(:)
         integer, allocatable :: old_to_new(:)
         integer :: n_used, j

         allocate(used(n_total), old_to_new(n_total))
         used = .false.
         do i = 1, iconf_data
            used(xyz_config_map(i)) = .true.
         end do
         n_used = count(used)
         allocate(xyz_configs_compact(n_used))
         old_to_new = 0
         j = 0
         do i = 1, n_total
            if (used(i)) then
               j = j + 1
               xyz_configs_compact(j) = xyz_configs(i)
               old_to_new(i) = j
            end if
         end do
         ! Update the map to new indices
         do i = 1, iconf_data
            xyz_config_map(i) = old_to_new(xyz_config_map(i))
         end do
         ! Replace the global array
         deallocate(xyz_configs)
         allocate(xyz_configs(n_used))
         xyz_configs(1:n_used) = xyz_configs_compact(1:n_used)
         xyz_n_configs_total = n_used
         deallocate(xyz_configs_compact, used, old_to_new)
         if (rangml == 0) then
            write(6, '("ML: XYZ configs compacted: ", i0, " of ", i0, " kept")') n_used, n_total
         end if
      end block

      if (rangml == 0) then
         write(6, '("ML: XYZ database prepared: ", i0, " total, ", i0, " train, ", i0, " test")') &
            iconf_data, iconf_data_train, iconf_data_test
      end if

      deallocate(class_list)

      _MLD_END_
   end subroutine xyz_prepare_database


   !---------------------------------------------------------------------------
   !> Read a single configuration from the XYZ file.
   !> This replaces read_poscar_sasha for XYZ mode.
   !> Only subrank==0 reads the file; data is broadcast to subworld.
   !---------------------------------------------------------------------------
   subroutine read_xyz_config(ifile)
      use ml_in_ndm_module, only: rangml, debug, im, imm, &
         fix_weighted_for_element, fix_weighted_for_element_3ch, &
         fix_weighted_for_element_ini, fix_weighted_for_element_3ch_ini, &
         weighted, weighted_3ch, &
         linvisible, fix_ch_elements_invisible, fix_no_of_elements_invisible
      use module_chemical_species, only: periodic_table_element, size_periodic_table, &
         fix_ch_elements, fix_no_of_elements
      use derived_types, only: config_real
      use module_db_poscar, only: fix_ref_energy_per_element
      use module_db_setup, only: db_path
      use math
      use my_mpi_subroutines, only: my_broadcast_char, my_broadcast_val_int, &
         my_broadcast_val_real, my_broadcast_vect_real, my_broadcast_vect_int
      use mld_unit
      use mld_logger
      use mld_string
      use mld_mpi
      use mld_subworld, only: subworld, subrank
      implicit none

      integer, intent(in) :: ifile

      integer :: ixyz  ! index into xyz_configs
      integer :: inp, ios, natoms, i, ij, ii, i_p, icnt
      integer :: im_local, nb_elements
      character(len=500) :: atom_line
      character(len=20) :: ctmp
      character(len=MAX_INFO_LEN) :: info_line
      character(len=2) :: ch_buffer
      real(kind_double), allocatable :: vect_buffer(:)
      real(kind_double) :: box(3,3), box_inv(3,3), st(6), volume, tmp_ene
      real(kind_double) :: E_total, E_free
      logical :: has_energy, has_forces, has_stress, has_lattice
      integer :: target_line, current_line
      ! For species parsing
      character(len=2), allocatable :: atom_species(:)
      real(kind_double), allocatable :: atom_pos(:,:), atom_force(:,:)
      ! For unique species
      character(len=2), allocatable :: element_list(:)
      integer, allocatable :: nspecies(:), ityp(:)
      real(kind_double), allocatable :: xp(:,:), xc(:,:), fp(:,:), l_spin(:,:)
      integer :: n_unique, found

      _NAMECURRENT_("read_xyz_config")
      _MLD_BEGIN_

      ! Map config_real index to xyz_configs index
      ixyz = xyz_config_map(ifile)

      natoms      = xyz_configs(ixyz)%natoms
      target_line = xyz_configs(ixyz)%file_line
      has_energy  = xyz_configs(ixyz)%has_energy
      has_forces  = xyz_configs(ixyz)%has_forces
      has_stress  = xyz_configs(ixyz)%has_stress
      has_lattice = xyz_configs(ixyz)%has_lattice
      E_total     = xyz_configs(ixyz)%energy
      E_free      = xyz_configs(ixyz)%free_energy
      box         = xyz_configs(ixyz)%lattice

      ! Remap stress from 3x3 row-major (xx,xy,xz,yx,yy,yz,zx,zy,zz) to Voigt (xx,yy,zz,xy,xz,yz)
      st(1) = xyz_configs(ixyz)%stress(1) ! xx
      st(2) = xyz_configs(ixyz)%stress(5) ! yy
      st(3) = xyz_configs(ixyz)%stress(9) ! zz
      st(4) = xyz_configs(ixyz)%stress(2) ! xy
      st(5) = xyz_configs(ixyz)%stress(3) ! xz
      st(6) = xyz_configs(ixyz)%stress(6) ! yz

      ! ---- Read atom data (subrank 0 only) ----
      allocate(atom_species(natoms))
      allocate(atom_pos(3,natoms))
      allocate(atom_force(3,natoms))
      atom_species(:) = '  '
      atom_pos = 0.d0
      atom_force = 0.d0

      if (subrank == 0) then
         open(newunit=inp, file=trim(db_path), status='old', action='read', access='stream', form='formatted', iostat=ios)
         if (ios /= 0) then
            write(6,*) 'ML: FATAL: cannot open XYZ file: ', trim(db_path)
            stop 'FATAL: cannot open XYZ file in read_xyz_config'
         end if

         ! Seek directly to first atom line (O(1) via stream byte offset)
         read(inp, '(a)', pos=xyz_configs(ixyz)%byte_offset, iostat=ios) atom_line
         if (ios /= 0) then
            write(6,*) 'ML: FATAL: unexpected EOF reading atoms for config ', ifile
            stop 'FATAL: unexpected EOF reading atoms in read_xyz_config'
         end if
         if (has_forces) then
            read(atom_line, *, iostat=ios) atom_species(1), &
               atom_pos(1,1), atom_pos(2,1), atom_pos(3,1), &
               atom_force(1,1), atom_force(2,1), atom_force(3,1)
         else
            read(atom_line, *, iostat=ios) atom_species(1), &
               atom_pos(1,1), atom_pos(2,1), atom_pos(3,1)
         end if
         if (ios /= 0) then
            write(6,*) 'ML: WARNING: parse error for atom 1 in config ', ifile
            atom_pos(:,1) = 0.d0 ; atom_force(:,1) = 0.d0
         end if

         ! Read remaining atom lines sequentially
         do i = 2, natoms
            read(inp, '(a)', iostat=ios) atom_line
            if (ios /= 0) then
               write(6,*) 'ML: FATAL: unexpected EOF reading atoms for config ', ifile
               stop 'FATAL: unexpected EOF reading atoms in read_xyz_config'
            end if
            if (has_forces) then
               read(atom_line, *, iostat=ios) atom_species(i), &
                  atom_pos(1,i), atom_pos(2,i), atom_pos(3,i), &
                  atom_force(1,i), atom_force(2,i), atom_force(3,i)
            else
               read(atom_line, *, iostat=ios) atom_species(i), &
                  atom_pos(1,i), atom_pos(2,i), atom_pos(3,i)
            end if
            if (ios /= 0) then
               write(6,*) 'ML: WARNING: parse error for atom ', i, ' in config ', ifile
               atom_pos(:,i) = 0.d0 ; atom_force(:,i) = 0.d0
            end if
         end do
         close(inp)
      end if

      ! ---- Broadcast atom data to subworld ----
      ! Species (broadcast each 2-char element)
      do i = 1, natoms
         ch_buffer = atom_species(i)
         call my_broadcast_char(ch_buffer, 0, subworld, codeml)
         atom_species(i) = ch_buffer
      end do

      ! Positions and forces (broadcast as vectors per component)
      allocate(vect_buffer(natoms))
      do ii = 1, 3
         vect_buffer(1:natoms) = atom_pos(ii,1:natoms)
         call my_broadcast_vect_real(vect_buffer, 0, subworld, codeml)
         atom_pos(ii,1:natoms) = vect_buffer(1:natoms)
      end do
      do ii = 1, 3
         vect_buffer(1:natoms) = atom_force(ii,1:natoms)
         call my_broadcast_vect_real(vect_buffer, 0, subworld, codeml)
         atom_force(ii,1:natoms) = vect_buffer(1:natoms)
      end do
      deallocate(vect_buffer)

      ! Broadcast box (9 values) and stress (6 values)
      allocate(vect_buffer(9))
      vect_buffer(1:3) = box(:,1)
      vect_buffer(4:6) = box(:,2)
      vect_buffer(7:9) = box(:,3)
      call my_broadcast_vect_real(vect_buffer, 0, subworld, codeml)
      box(:,1) = vect_buffer(1:3)
      box(:,2) = vect_buffer(4:6)
      box(:,3) = vect_buffer(7:9)
      deallocate(vect_buffer)

      allocate(vect_buffer(6))
      vect_buffer(1:6) = st(1:6)
      call my_broadcast_vect_real(vect_buffer, 0, subworld, codeml)
      st(1:6) = vect_buffer(1:6)
      deallocate(vect_buffer)

      call my_broadcast_val_real(E_total, 0, subworld, codeml)
      call my_broadcast_val_real(E_free,  0, subworld, codeml)

      ! ---- Determine unique species (sorted alphabetically like convert_db.py) ----
      n_unique = 0
      allocate(element_list(natoms))
      element_list(:) = '  '
      do i = 1, natoms
         found = 0
         do ij = 1, n_unique
            if (element_list(ij) == atom_species(i)) then
               found = ij
               exit
            end if
         end do
         if (found == 0) then
            n_unique = n_unique + 1
            element_list(n_unique) = atom_species(i)
         end if
      end do
      ! Sort species alphabetically (simple bubble sort, n_unique is small)
      call sort_species(element_list, n_unique)

      nb_elements = n_unique
      im_local = natoms

      ! Count atoms per species
      allocate(nspecies(nb_elements))
      nspecies = 0
      do i = 1, natoms
         do ij = 1, nb_elements
            if (atom_species(i) == element_list(ij)) then
               nspecies(ij) = nspecies(ij) + 1
               exit
            end if
         end do
      end do

      ! Build type array and reorder atoms by species (POSCAR convention)
      allocate(ityp(im_local))
      allocate(xp(3,im_local), xc(3,im_local), fp(3,im_local), l_spin(3,im_local))
      l_spin = 0.d0
      icnt = 0
      do ij = 1, nb_elements
         do i = 1, natoms
            if (atom_species(i) == element_list(ij)) then
               icnt = icnt + 1
               ityp(icnt) = ij
               xp(:,icnt) = atom_pos(:,i)
               fp(:,icnt) = atom_force(:,i)
            end if
         end do
      end do

      ! Positions are Cartesian in XYZ; compute crystal coords.
      ! has_lattice alone is the correct test: it's set purely from whether
      ! Lattice="..." was present and parsed as 9 numbers (see
      ! parse_xyz_info_line), unlike a diagonal-magnitude check, which
      ! misfires on valid zero-diagonal cells (e.g. the standard FCC
      ! primitive-cell convention a=(0,L,L), b=(L,0,L), c=(L,L,0)).
      if (has_lattice) then
         call matinv_gen(box, box_inv)
         xc(1:3, 1:im_local) = matmul(box_inv(1:3,1:3), xp(1:3,1:im_local))
      else
         ! No lattice (isolated molecule): set box to large cube
         box = 0.d0
         box(1,1) = 100.d0
         box(2,2) = 100.d0
         box(3,3) = 100.d0
         call matinv_gen(box, box_inv)
         xc(1:3, 1:im_local) = matmul(box_inv(1:3,1:3), xp(1:3,1:im_local))
      end if

      ! ---- Update global im/imm ----
      if (im_local > imm) then
         imm = im_local
         if (debug) write(6,'("ML: WARNING imm increased to ",i7," in read_xyz_config")') imm
      end if
      im = im_local

      ! ---- Fill config_real arrays (same pattern as read_poscar_sasha) ----
      if (allocated(config_real(ifile)%ref_energy_per_element)) deallocate(config_real(ifile)%ref_energy_per_element)
      allocate(config_real(ifile)%ref_energy_per_element(nb_elements))
      if (allocated(config_real(ifile)%mass_per_type)) deallocate(config_real(ifile)%mass_per_type)
      allocate(config_real(ifile)%mass_per_type(nb_elements))
      if (allocated(config_real(ifile)%Z_per_type)) deallocate(config_real(ifile)%Z_per_type)
      allocate(config_real(ifile)%Z_per_type(nb_elements))
      if (allocated(config_real(ifile)%covalent_radius_per_type)) deallocate(config_real(ifile)%covalent_radius_per_type)
      allocate(config_real(ifile)%covalent_radius_per_type(nb_elements))
      if (allocated(config_real(ifile)%fix_type_poscar_to_periodic)) deallocate(config_real(ifile)%fix_type_poscar_to_periodic)
      allocate(config_real(ifile)%fix_type_poscar_to_periodic(nb_elements))
      if (allocated(config_real(ifile)%itype_to_global)) deallocate(config_real(ifile)%itype_to_global)
      allocate(config_real(ifile)%itype_to_global(nb_elements))

      ! Assign mass/Z/radius from periodic table
      do ij = 1, nb_elements
         do i_p = 1, size_periodic_table
            if (element_list(ij) == periodic_table_element(i_p)%symbol) then
               config_real(ifile)%fix_type_poscar_to_periodic(ij) = periodic_table_element(i_p)%Z
               config_real(ifile)%mass_per_type(ij) = periodic_table_element(i_p)%mass
               config_real(ifile)%Z_per_type(ij) = periodic_table_element(i_p)%Z
               config_real(ifile)%covalent_radius_per_type(ij) = periodic_table_element(i_p)%covalent_radius
               exit
            end if
         end do
      end do

      ! Assign reference energies and type mapping
      if (weighted) then
         do ij = 1, nb_elements
            icnt = 0
            do i_p = 1, fix_no_of_elements
               if (fix_ch_elements(i_p) == element_list(ij)) then
                  config_real(ifile)%itype_to_global(ij) = i_p
                  icnt = icnt + 1
                  config_real(ifile)%ref_energy_per_element(ij) = fix_ref_energy_per_element(i_p)
               end if
            end do
            if (icnt == 0) then
               write(6,*) "ML: Unknown element in XYZ database: ", element_list(ij)
               write(6,*) "ML: Please update chemical_elements in input file."
               stop 'fatal in read_xyz_config, unknown element (weighted)'
            end if
         end do
         if (linvisible) then
            if (allocated(config_real(ifile)%invisible_per_type)) deallocate(config_real(ifile)%invisible_per_type)
            allocate(config_real(ifile)%invisible_per_type(nb_elements))
            config_real(ifile)%invisible_per_type(:) = .false.
            do ij = 1, nb_elements
               do i_p = 1, fix_no_of_elements_invisible
                  if (fix_ch_elements_invisible(i_p) == element_list(ij)) then
                     config_real(ifile)%invisible_per_type(ij) = .true.
                  end if
               end do
            end do
         end if
      else
         do ij = 1, nb_elements
            icnt = 0
            do i_p = 1, fix_no_of_elements
               if (fix_ch_elements(i_p) == element_list(ij)) then
                  icnt = icnt + 1
                  config_real(ifile)%ref_energy_per_element(ij) = fix_ref_energy_per_element(i_p)
                  config_real(ifile)%itype_to_global(ij) = i_p
               end if
            end do
            if (icnt == 0) then
               write(6,*) "ML: Unknown element in XYZ database: ", element_list(ij)
               write(6,*) "ML: Please update chemical_elements."
               stop 'fatal in read_xyz_config, unknown element'
            end if
         end do
      end if

      if (weighted) then
         fix_weighted_for_element_ini(:) = fix_weighted_for_element(:)
         if (weighted_3ch) then
            fix_weighted_for_element_3ch_ini(:) = fix_weighted_for_element_3ch(:)
         end if
      end if

      ! ---- Store into config_real ----
      if (allocated(config_real(ifile)%itype)) deallocate(config_real(ifile)%itype)
      allocate(config_real(ifile)%itype(im_local))
      if (allocated(config_real(ifile)%itype_db)) deallocate(config_real(ifile)%itype_db)
      allocate(config_real(ifile)%itype_db(im_local))
      if (allocated(config_real(ifile)%pos_cart)) deallocate(config_real(ifile)%pos_cart)
      allocate(config_real(ifile)%pos_cart(3,im_local))
      if (allocated(config_real(ifile)%pos_crst)) deallocate(config_real(ifile)%pos_crst)
      allocate(config_real(ifile)%pos_crst(3,im_local))
      if (allocated(config_real(ifile)%force)) deallocate(config_real(ifile)%force)
      allocate(config_real(ifile)%force(3,im_local))
      if (allocated(config_real(ifile)%atomic_spin)) deallocate(config_real(ifile)%atomic_spin)
      allocate(config_real(ifile)%atomic_spin(3,im_local))

      config_real(ifile)%ntypes = nb_elements
      config_real(ifile)%nat = im_local
      config_real(ifile)%itype(1:im_local) = ityp(1:im_local)

      do ii = 1, im_local
         config_real(ifile)%itype_db(ii) = config_real(ifile)%itype_to_global(ityp(ii))
      end do

      config_real(ifile)%pos_crst(1:3,1:im_local) = xc(1:3,1:im_local)
      config_real(ifile)%pos_cart(1:3,1:im_local)  = xp(1:3,1:im_local)
      config_real(ifile)%force(1:3,1:im_local)     = fp(1:3,1:im_local)

      ! Subtract reference energy
      tmp_ene = 0.d0
      do ii = 1, im_local
         tmp_ene = tmp_ene + config_real(ifile)%ref_energy_per_element(ityp(ii))
      end do
      config_real(ifile)%ref_energy = tmp_ene
      config_real(ifile)%energy(1) = E_total - tmp_ene
      config_real(ifile)%energy(2) = E_free  - tmp_ene
      config_real(ifile)%energy(3) = E_free  - tmp_ene

      config_real(ifile)%cell = box
      config_real(ifile)%stress(1:6) = st(1:6)
      config_real(ifile)%spin = 0
      config_real(ifile)%atomic_spin = l_spin

      config_real(ifile)%has_energy = has_energy .and. config_real(ifile)%has_energy
      config_real(ifile)%has_force  = has_forces .and. config_real(ifile)%has_force
      config_real(ifile)%has_stress = has_stress .and. config_real(ifile)%has_stress

      call calc_volume(config_real(ifile)%cell(:,1), config_real(ifile)%cell(:,2), &
                        config_real(ifile)%cell(:,3), volume)
      config_real(ifile)%volume = volume

      ! Cleanup
      deallocate(atom_species, atom_pos, atom_force)
      deallocate(element_list, nspecies, ityp, xp, xc, fp, l_spin)

      _MLD_END_
   end subroutine read_xyz_config


   !---------------------------------------------------------------------------
   ! PRIVATE helper routines
   !---------------------------------------------------------------------------

   !> Parse the info line of an extended XYZ config
   subroutine parse_xyz_info_line(line, cfg)
      implicit none
      character(len=*), intent(in) :: line
      type(xyz_config_info), intent(inout) :: cfg
      integer :: i, j, ios
      character(len=MAX_INFO_LEN) :: val_str
      real(kind_double) :: sv(9)
      logical :: found_ene

      ! Defaults
      cfg%has_energy  = .false.
      cfg%has_forces  = .false.
      cfg%has_stress  = .false.
      cfg%has_lattice = .false.
      cfg%energy      = 0.d0
      cfg%free_energy = 0.d0
      cfg%stress      = 0.d0
      cfg%lattice     = 0.d0
      cfg%class_id    = '01'

      ! ---- Check Properties for forces ----
      i = index(line, 'Properties=')
      if (i > 0) then
         val_str = line(i+11:)
         if (index(val_str, 'forces') > 0 .or. index(val_str, 'Forces') > 0) then
            cfg%has_forces = .true.
         end if
      end if

      ! ---- Parse Lattice="..." ----
      i = index(line, 'Lattice="')
      if (i > 0) then
         j = index(line(i+9:), '"')
         if (j > 0) then
            val_str = line(i+9 : i+9+j-2)
            read(val_str, *, iostat=ios) cfg%lattice(1,1), cfg%lattice(2,1), cfg%lattice(3,1), &
                                          cfg%lattice(1,2), cfg%lattice(2,2), cfg%lattice(3,2), &
                                          cfg%lattice(1,3), cfg%lattice(2,3), cfg%lattice(3,3)
            if (ios == 0) cfg%has_lattice = .true.
         end if
      end if

      ! ---- Parse the training energy target ----
      ! Prefer atomization_energy= (energy relative to isolated-atom
      ! references - the physically meaningful ML fitting target) when
      ! present, e.g. the mad-1.5 dataset. Fall back to the raw energy= key
      ! (careful: must not match free_energy=) for datasets that only have
      ! that, e.g. HEA25.extxyz.
      call extract_kv_real(line, 'atomization_energy=', cfg%energy, found_ene)
      if (.not. found_ene) then
         call extract_kv_real(line, 'energy=', cfg%energy, found_ene)
      end if
      cfg%has_energy = found_ene

      ! ---- Parse free_energy=... ----
      call extract_kv_real(line, 'free_energy=', cfg%free_energy, found_ene)
      if (.not. found_ene) cfg%free_energy = cfg%energy

      ! ---- Parse stress="..." ----
      i = index(line, 'stress="')
      if (i > 0) then
         j = index(line(i+8:), '"')
         if (j > 0) then
            val_str = line(i+8 : i+8+j-2)
            read(val_str, *, iostat=ios) sv(1:9)
            if (ios == 0) then
               cfg%stress(1:9) = sv(1:9)
               cfg%has_stress = .true.
            end if
         end if
      end if

      ! ---- Parse class=... ----
      call extract_kv_string(line, 'class=', val_str, ios)
      if (ios == 0) then
         val_str = adjustl(val_str)
         j = len_trim(val_str)
         if (j == 1) then
            cfg%class_id = '0'//val_str(1:1)
         else if (j >= 2) then
            cfg%class_id = val_str(1:2)
         end if
      else
         ! No class= key: fall back to subset=, auto-registering each
         ! distinct name as a new class in first-appearance order. This is
         ! completely dataset-independent - no dataset-specific table.
         call extract_kv_string(line, 'subset=', val_str, ios)
         if (ios == 0) then
            cfg%class_id = auto_register_subset_class(trim(adjustl(val_str)))
         end if
      end if

   end subroutine parse_xyz_info_line


   !> Auto-register a subset= name, assigning class ids in first-appearance
   !> order starting at '01'. Returns the same id for a name seen again
   !> later in the file. Dataset-independent: no fixed name table.
   function auto_register_subset_class(subset_name) result(class_id)
      implicit none
      character(len=*), intent(in) :: subset_name
      character(len=2) :: class_id
      integer :: k

      do k = 1, n_auto_subset
         if (trim(auto_subset_name(k)) == trim(subset_name)) then
            class_id = auto_subset_class(k)
            return
         end if
      end do

      n_auto_subset = n_auto_subset + 1
      if (n_auto_subset > size(auto_subset_name)) then
         write(6,*) 'ML: FATAL: more than ', size(auto_subset_name), &
            ' distinct auto-registered subset= names in XYZ file'
         stop 'too many auto-registered subset classes in XYZ'
      end if
      auto_subset_name(n_auto_subset) = subset_name
      write(auto_subset_class(n_auto_subset), '(i2.2)') n_auto_subset
      class_id = auto_subset_class(n_auto_subset)
   end function auto_register_subset_class


   !> Extract a real-valued key from "key=value" in a line.
   !> Ensures word boundary so "energy=" doesn't match "free_energy=".
   subroutine extract_kv_real(line, key, val, found)
      implicit none
      character(len=*), intent(in) :: line, key
      real(kind_double), intent(out) :: val
      logical, intent(out) :: found
      integer :: i, j, ios, l, start_pos, klen

      found = .false.
      val = 0.d0
      l = len_trim(line)
      klen = len_trim(key)
      start_pos = 1

      ! Search for a match at a proper word boundary; a match that fails the
      ! boundary check (e.g. "energy=" inside "atomization_energy=") does not
      ! stop the search - keep looking for a later, real occurrence of key.
      do
         if (start_pos > l) return
         i = index(line(start_pos:l), trim(key))
         if (i == 0) return
         i = i + start_pos - 1

         if (i == 1 .or. line(i-1:i-1) == ' ' .or. line(i-1:i-1) == '"' .or. line(i-1:i-1) == ',') then
            j = i + klen
            if (j > l) return
            read(line(j:), *, iostat=ios) val
            if (ios == 0) found = .true.
            return
         end if

         start_pos = i + klen
      end do
   end subroutine extract_kv_real


   !> Extract a string-valued key from "key=value" in a line (no quotes)
   subroutine extract_kv_string(line, key, val, ios)
      implicit none
      character(len=*), intent(in) :: line, key
      character(len=*), intent(out) :: val
      integer, intent(out) :: ios
      integer :: i, j, k, l

      ios = 1  ! not found
      val = ' '
      l = len_trim(line)
      i = index(line, trim(key))
      if (i == 0) return

      j = i + len_trim(key)
      if (j > l) return

      ! Find end of value (next space)
      k = index(line(j:), ' ')
      if (k == 0) then
         val = line(j:l)
      else
         val = line(j:j+k-2)
      end if
      ios = 0
   end subroutine extract_kv_string


   !> Find a class in the list or add it; return its index
   function find_or_add_class(class_id, labels, counts, n) result(idx)
      implicit none
      character(len=2), intent(in) :: class_id
      character(len=2), intent(inout) :: labels(100)
      integer, intent(inout) :: counts(100)
      integer, intent(inout) :: n
      integer :: idx, i

      do i = 1, n
         if (labels(i) == class_id) then
            idx = i
            return
         end if
      end do
      ! Not found, add new
      n = n + 1
      if (n > 100) then
         write(6,*) 'ML: FATAL: more than 100 classes in XYZ file'
         stop 'too many classes in XYZ'
      end if
      labels(n) = class_id
      counts(n) = 0
      idx = n
   end function find_or_add_class


   !> Sort species names alphabetically (bubble sort, n is small)
   subroutine sort_species(arr, n)
      implicit none
      integer, intent(in) :: n
      character(len=2), intent(inout) :: arr(n)
      character(len=2) :: tmp
      integer :: i, j
      logical :: swapped

      do i = 1, n-1
         swapped = .false.
         do j = 1, n-i
            if (lgt(arr(j), arr(j+1))) then
               tmp = arr(j)
               arr(j) = arr(j+1)
               arr(j+1) = tmp
               swapped = .true.
            end if
         end do
         if (.not. swapped) exit
      end do
   end subroutine sort_species

end module module_extxyz


!---------------------------------------------------------------------------
!> Check whether any pair of atoms in a configuration is closer than rmin.
!> Shared by both the POSCAR and the XYZ database-reading paths so that
!> drop_short_dist behaves identically regardless of the input format.
!> Uses a single nearest-image convention (fractional round-to-nearest),
!> matching the reference Python filtering scripts used to build the
!> "noshort" training databases.
!---------------------------------------------------------------------------
subroutine config_has_short_pair(nat, pos_cart, cell, rmin, too_close)
   use module_kind_variables, only: kind_double
   use math, only: matinv_gen
   implicit none

   integer, intent(in)  :: nat
   real(kind_double), intent(in)  :: pos_cart(3, nat)
   real(kind_double), intent(in)  :: cell(3, 3)
   real(kind_double), intent(in)  :: rmin
   logical, intent(out) :: too_close

   real(kind_double) :: cell_inv(3, 3)
   real(kind_double) :: frac(3, nat)
   real(kind_double) :: dfrac(3), dcart(3), r2, rmin2
   integer :: ia, ja

   too_close = .false.
   if (nat <= 1 .or. rmin <= 0.d0) return
   rmin2 = rmin*rmin

   call matinv_gen(cell, cell_inv)
   do ia = 1, nat
      frac(:, ia) = matmul(cell_inv, pos_cart(:, ia))
   end do

   do ia = 1, nat - 1
      do ja = ia + 1, nat
         dfrac = frac(:, ia) - frac(:, ja)
         dfrac = dfrac - nint(dfrac)
         dcart = matmul(cell, dfrac)
         r2 = dot_product(dcart, dcart)
         if (r2 < rmin2) then
            too_close = .true.
            return
         end if
      end do
   end do
end subroutine config_has_short_pair


! ============================================================================
! poscar_scan_pos_cell has an allocatable dummy argument, which requires an
! explicit interface at every call site; wrapping it in a module (instead of
! leaving it external like the other standalone subroutines in this file)
! gives the compiler that interface automatically.
! ============================================================================
module module_poscar_short_scan
   implicit none
   private
   public :: poscar_scan_pos_cell
contains

   !---------------------------------------------------------------------------
   !> Minimal, rank-0-only reader of just the lattice + Cartesian positions
   !> from a Milady POSCAR file, used by the drop_short_dist pre-scan so that
   !> filtering can happen before random selection. No MPI, no element/weight
   !> bookkeeping - read_poscar_sasha remains the authoritative full reader
   !> used once a config has actually been selected.
   !---------------------------------------------------------------------------
subroutine poscar_scan_pos_cell(name_file, ok, nat, pos_cart, cell)
   use module_kind_variables, only: kind_double
   implicit none
   character(len=*), intent(in)  :: name_file
   logical, intent(out) :: ok
   integer, intent(out) :: nat
   real(kind_double), allocatable, intent(out) :: pos_cart(:, :)
   real(kind_double), intent(out) :: cell(3, 3)

   integer :: inp, ios, nb_elements, nitype, ij, i, icnt, i_pos_form
   integer :: nitems2, itest_there_is_a_number
   real(kind_double) :: alat, e1, e2, e3
   character(len=3) :: EFS_tag
   character(len=500) :: line
   character(len=80) :: dummy, trimdummy
   character(len=2), allocatable :: elem(:)
   integer, allocatable :: mass_dummy(:), nspecies(:)
   logical :: lexist

   ok = .false.
   nat = 0
   cell = 0.d0

   inquire (file=trim(name_file), exist=lexist)
   if (.not. lexist) return
   open (newunit=inp, file=trim(name_file), status='old', action='read', iostat=ios)
   if (ios /= 0) return

   read (inp, *, iostat=ios) EFS_tag, nb_elements
   if (ios /= 0 .or. nb_elements < 1) then; close (inp); return; end if
   backspace (inp)

   allocate (elem(nb_elements), mass_dummy(nb_elements))
   read (inp, *, iostat=ios) EFS_tag, nb_elements, (elem(ij), mass_dummy(ij), ij=1, nb_elements), e1, e2, e3
   if (ios /= 0) then; close (inp); return; end if

   read (inp, *, iostat=ios) alat
   if (ios /= 0) then; close (inp); return; end if
   read (inp, *, iostat=ios) cell(1:3, 1)
   read (inp, *, iostat=ios) cell(1:3, 2)
   read (inp, *, iostat=ios) cell(1:3, 3)
   if (ios /= 0) then; close (inp); return; end if
   cell = alat*cell

   read (inp, '(a)', iostat=ios) line
   if (ios /= 0) then; close (inp); return; end if
   if (itest_there_is_a_number(line) == 0) then
      nitype = nitems2(line)   ! VASP5 species-symbol line, already consumed
   else
      nitype = nitems2(line)
      backspace (inp)          ! VASP4: this was actually the counts line
   end if
   if (nitype /= nb_elements) then; close (inp); return; end if

   allocate (nspecies(nitype))
   read (inp, *, iostat=ios) (nspecies(i), i=1, nitype)
   if (ios /= 0) then; close (inp); return; end if
   read (inp, '(a)', iostat=ios) dummy
   if (ios /= 0) then; close (inp); return; end if
   trimdummy = trim(adjustl(dummy))
   if (trimdummy(1:1) == 'd' .or. trimdummy(1:1) == 'D') then
      i_pos_form = 1
   else if (trimdummy(1:1) == 'c' .or. trimdummy(1:1) == 'C') then
      i_pos_form = 0
   else
      close (inp); return
   end if

   nat = sum(nspecies)
   allocate (pos_cart(3, nat))
   icnt = 0
   do ij = 1, nitype
      do i = 1, nspecies(ij)
         icnt = icnt + 1
         read (inp, *, iostat=ios) pos_cart(1:3, icnt)
         if (ios /= 0) then; close (inp); return; end if
      end do
   end do
   close (inp)

   if (i_pos_form == 1) pos_cart = matmul(cell, pos_cart)   ! Direct -> Cartesian

   ok = .true.
end subroutine poscar_scan_pos_cell

end module module_poscar_short_scan


!---------------------------------------------------------------------------
!> Apply drop_short_dist filtering to the POSCAR database path: for every
!> class in db_model, pre-scan all candidate POSCAR files (rank 0 only,
!> broadcast to the rest), drop the ones with an interatomic distance below
!> the threshold, and compact no_total/no_start/no_selec so the subsequent
!> random selection (rks2, in prepare_name_file_from_db) only ever sees
!> surviving configs - mirroring how the XYZ path already excludes dropped
!> configs before selection. No-op when drop_short_dist <= 0.
!---------------------------------------------------------------------------
subroutine filter_short_distance_poscar()
   use module_kind_variables, only: kind_double
   use module_db_setup, only: db_path, drop_short_dist, db_train, db_test, &
      iconf_data, iconf_data_train, iconf_data_test
   use derived_types, only: db_model, config_real, config_desc
   use module_poscar_short_scan, only: poscar_scan_pos_cell
   use mld_mpi, only: mpi_comm_mld, rangml, mld_ierror
   use mpi
   implicit none

   integer :: i, k, n_candidates, n_kept, n_dropped, n_unreadable, iconf_temp
   integer, allocatable :: kept_tmp(:)
   character(len=180) :: ftemp
   character(len=7) :: cnumber_str
   integer, parameter :: number0 = 1000000
   logical :: ok, too_close
   integer :: nat
   real(kind_double), allocatable :: pos_cart(:, :)
   real(kind_double) :: cell(3, 3)

   if (drop_short_dist <= 0.d0) return

   do i = 1, size(db_model)
      n_candidates = db_model(i)%no_total - db_model(i)%no_start + 1
      if (allocated(kept_tmp)) deallocate (kept_tmp)
      allocate (kept_tmp(n_candidates))
      kept_tmp = 0
      n_kept = 0
      n_dropped = 0
      n_unreadable = 0

      if (rangml == 0) then
         do k = db_model(i)%no_start, db_model(i)%no_total
            write (cnumber_str, '(i7)') number0 + k
            ftemp = trim(adjustl(db_path))//db_model(i)%class//"_"//db_model(i)%klm//"_"//cnumber_str(2:7)//'.poscar'
            call poscar_scan_pos_cell(trim(ftemp), ok, nat, pos_cart, cell)
            if (.not. ok) then
               n_unreadable = n_unreadable + 1
               cycle
            end if
            call config_has_short_pair(nat, pos_cart, cell, drop_short_dist, too_close)
            if (too_close) then
               n_dropped = n_dropped + 1
               cycle
            end if
            n_kept = n_kept + 1
            kept_tmp(n_kept) = k
         end do
         if (n_unreadable > 0) then
            write (6, '("ML: WARNING: class ", a, " klm ", a, ": ", i0, &
               &" POSCAR file(s) could not be read during drop_short_dist pre-scan, excluded from the candidate pool")') &
               db_model(i)%class, db_model(i)%klm, n_unreadable
         end if
         write (6, '("ML: drop_short_dist = ", f10.5, " Angstrom: class ", a, " klm ", a, &
            &": kept ", i0, " of ", i0, " (dropped ", i0, ")")') &
            drop_short_dist, db_model(i)%class, db_model(i)%klm, n_kept, n_candidates, n_dropped
      end if

      call MPI_BCAST(n_kept, 1, MPI_INTEGER, 0, mpi_comm_mld, mld_ierror)
      call MPI_BCAST(kept_tmp, n_candidates, MPI_INTEGER, 0, mpi_comm_mld, mld_ierror)

      if (allocated(db_model(i)%kept_file_number)) deallocate (db_model(i)%kept_file_number)
      allocate (db_model(i)%kept_file_number(max(n_kept, 1)))
      db_model(i)%kept_file_number(1:n_kept) = kept_tmp(1:n_kept)

      db_model(i)%no_total = n_kept
      db_model(i)%no_start = 1
      if (db_model(i)%no_selec > db_model(i)%no_total) then
         if (rangml == 0) write (6, '("ML: WARNING: class ", a, " klm ", a, ": no_selec=", i0, &
            &" exceeds filtered no_total=", i0, ", clamping")') &
            db_model(i)%class, db_model(i)%klm, db_model(i)%no_selec, db_model(i)%no_total
         db_model(i)%no_selec = db_model(i)%no_total
      end if
   end do
   if (allocated(kept_tmp)) deallocate (kept_tmp)

   ! Recompute iconf_data/train/test after filtering and re-allocate
   ! config_real/config_desc/db_train/db_test to the corrected sizes
   ! (read_db_file sized them from the pre-filter no_total/no_selec).
   iconf_temp = 0
   do i = 1, size(db_model)
      iconf_temp = iconf_temp + db_model(i)%no_selec
   end do
   iconf_data_train = iconf_temp

   iconf_temp = 0
   do i = 1, size(db_model)
      iconf_temp = iconf_temp + db_model(i)%no_total - db_model(i)%no_start + 1
   end do
   iconf_data = iconf_temp
   iconf_data_test = iconf_data - iconf_data_train

   if (allocated(config_real)) deallocate (config_real); allocate (config_real(iconf_data))
   if (allocated(config_desc)) deallocate (config_desc); allocate (config_desc(iconf_data))
   if (allocated(db_train)) deallocate (db_train); allocate (db_train(iconf_data_train))
   if (allocated(db_test)) deallocate (db_test); allocate (db_test(max(iconf_data_test, 0)))

   if (rangml == 0) then
      write (6, '("ML: drop_short_dist: database recomputed after filtering: ", i0, &
         &" total, ", i0, " train, ", i0, " test")') iconf_data, iconf_data_train, iconf_data_test
   end if
end subroutine filter_short_distance_poscar


! ============================================================================
! Standalone subroutines for database preparation and POSCAR reading
! ============================================================================

subroutine prepare_database()

   use module_db_setup, only: db_file
   use module_extxyz, only: db_xyz, xyz_scan_database, xyz_prepare_database
   use module_json, only: db_json, json_scan_database, json_prepare_database
   use mld_unit
   use mld_logger
   use mld_mpi

   implicit none

   _NAMECURRENT_("prepare_database")


   integer  :: n_db_class, iunit, ierr
   logical  :: ok, lexist

   _MLD_BEGIN_

   ! ---- Extended XYZ mode: scan the XYZ file and build database from it ----
   if (db_xyz) then
      call xyz_scan_database()
      call xyz_prepare_database()
      _MLD_END_
      return
   end if

   ! ---- MPtrj-style JSON mode: scan the JSON file and build database from it ----
   if (db_json) then
      call json_scan_database()
      call json_prepare_database()
      _MLD_END_
      return
   end if

   inquire(file=db_file, exist=lexist)
   if (.not. lexist) then
      call log_critical('prepare_database: db_file not found: '//trim(db_file))
      call log_critical('prepare_database: please provide the db_file')
      call mld_mpi_finalize('prepare_database: please provide the db_file')
      stop
   end if
   open (newunit=iunit, file=db_file, status='old', action='read', iostat=ierr)
   if (ierr /= 0) then
      call log_critical('prepare_database: error opening db_file: '//trim(db_file))
      call log_critical('prepare_database: check file permissions')
      call mld_mpi_finalize('prepare_database: check file permissions')
      stop
   end if
   
   ! reading how many lines are in db_file
   call read_n_db_class(iunit, n_db_class)
   ! allocate and fill the db_model object
   call read_db_file(iunit, n_db_class)
   ! drop configs with an interatomic distance below drop_short_dist, if
   ! requested; recomputes no_total/no_selec and re-sizes config_real/etc.
   call filter_short_distance_poscar()
   !if (ml_type == ml_type_analysis)
   call get_distinct_classes()
   close(iunit)
   call prepare_name_file_from_db()
   _MLD_END_
end subroutine prepare_database



subroutine get_distinct_classes()
   use ml_in_ndm_module, only: classes_full_for_sigma
   use module_error_by_class, only: number_of_distinct_classes
   use derived_types, only: db_model
   use mld_logger

   implicit none

   integer  :: k, i
   character(len=2), dimension(size(db_model, 1))     :: cvector, ctmp

   _NAMECURRENT_("get_distinct_classes")

   _MLD_BEGIN_
   ! initialize to avoid uninitialized trailing bytes
   cvector = '  '
   ctmp    = '  '
   do i = 1, size(db_model, 1)
      ctmp(i) = db_model(i)%class
   end do

   k = 1
   cvector(1) = ctmp(1)
   do i = 2, size(ctmp, 1)
      ! if the number already exist in res check next
      if (any(cvector == ctmp(i))) cycle
      ! No match found so add it to the output
      k = k + 1
      cvector(k) = ctmp(i)
   end do

   if (allocated(classes_full_for_sigma)) deallocate (classes_full_for_sigma); allocate (classes_full_for_sigma(k))
   classes_full_for_sigma(:) = cvector(1:k)

   number_of_distinct_classes = k

   _MLD_END_
end subroutine get_distinct_classes



subroutine read_n_db_class(inp, nlines)
   ! Read how many lines are in the inp file.
   ! The lines that have # in the firsrt character are ignored
   ! Input:
   !          inp (the file number)
   ! Output:
   !          nlines number of lines

   use ml_in_ndm_module, only: debug, rangml
   use mld_logger

   implicit none

   integer, intent(in)  :: inp
   integer, intent(out) :: nlines
   character(len=80)    :: ctmp
   integer  :: io

   _NAMECURRENT_("read_n_db_class")



   _MLD_BEGIN_
   nlines = 0
   DO
      read (inp, '(a)', iostat=io) ctmp
      !write (*,*) ctmp, ctmp(1:1), io
      if (io /= 0) exit
      if (ctmp(1:1) == '#') cycle
      nlines = nlines + 1
   END DO

   if (debug) then
      if (rangml == 0) write (6, '("ML: in read_n_db_class number of active lines in db file ", i5)') nlines
   end if

   rewind (inp)

   _MLD_END_
end subroutine read_n_db_class



subroutine read_db_file(inp, n_db_class)
   ! Read the n_db_class lines of the db file
   ! Input:
   !        inp (the file)
   !        n_db_class  (no of lines)
   ! Output:
   !        db_model object
   !        iconf_data, iconf_data_test, iconf_data_train
   !        allocate of config_real, config_desc, db_test, db_train

   use ml_in_ndm_module, only: rangml, debug
   use module_db_setup, only : selection_type, selection_type_first_start, db_train, db_test, &
      db_file,  iconf_data, iconf_data_train, iconf_data_test
   use module_optimization, only: optimize_weights_db
   use derived_types, only: db_model, config_real, config_desc
   use mld_logger
   use mld_mpi

   implicit none

   integer, intent(in)  :: inp, n_db_class
   character(len=80)    :: ctmp
   character(len=1)     :: c1, c2, c3
   integer  :: io, i, iconf_temp

   _NAMECURRENT_("read_db_file")



   _MLD_BEGIN_
   if (allocated(db_model)) deallocate (db_model); allocate (db_model(n_db_class))
   i = 0
   DO

      read (inp, '(a)', iostat=io) ctmp
      ! write (*,*) ctmp(1:1), ctmp
      if (io /= 0) exit
      if (ctmp(1:1) == '#') cycle
      backspace (inp)
      i = i + 1

      if (.not. (selection_type == selection_type_first_start)) then

         if (optimize_weights_db) then
            read (inp, *) db_model(i)%class, db_model(i)%klm, db_model(i)%no_total, db_model(i)%no_selec, c1, c2, c3, &
               db_model(i)%w_e, db_model(i)%w_f, db_model(i)%w_s, &
               db_model(i)%w_e_end, db_model(i)%w_f_end, db_model(i)%w_s_end
         else
            read (inp, *) db_model(i)%class, db_model(i)%klm, db_model(i)%no_total, db_model(i)%no_selec, c1, c2, c3, db_model(i)%w_e, db_model(i)%w_f, db_model(i)%w_s
         end if
         db_model(i)%no_start = 1

         if (db_model(i)%no_total .lt. db_model(i)%no_selec) then
            if (rangml == 0) then
               write (6, *) 'Some problems in read_db_file', db_file
               write (6, *) 'In the class', db_model(i)%class, 'and the KLM', db_model(i)%klm
               write (6, *) 'The number of total files', db_model(i)%no_total, 'is lower  than the no of selected files ', db_model(i)%no_selec
            end if
            call mld_mpi_finalize("no of files no_total and no_selec in read_db_file")
            stop
         end if

      else

         if (optimize_weights_db) then
            read (inp, *) db_model(i)%class, db_model(i)%klm, db_model(i)%no_total, db_model(i)%no_selec, db_model(i)%no_start, c1, c2, c3, &
               db_model(i)%w_e, db_model(i)%w_f, db_model(i)%w_s, &
               db_model(i)%w_e_end, db_model(i)%w_f_end, db_model(i)%w_s_end
         else
            read (inp, *) db_model(i)%class, db_model(i)%klm, db_model(i)%no_total, db_model(i)%no_selec, db_model(i)%no_start, c1, c2, c3, db_model(i)%w_e, db_model(i)%w_f, db_model(i)%w_s
         end if
         if (db_model(i)%no_total .lt. (db_model(i)%no_selec + db_model(i)%no_start - 1)) then
            if (rangml == 0) then
               write (6, *) 'Some problems in read_db_file', db_file
               write (6, *) 'In the class', db_model(i)%class, 'and the KLM', db_model(i)%klm
               write (6, *) 'The number of total files', db_model(i)%no_total, 'is lower than the no of selected files  + starting no of file ', db_model(i)%no_selec, db_model(i)%no_start
            end if
            call mld_mpi_finalize("no of files no_total and no_selec in read_db_file")
            stop
         end if

      end if                  ! selection_type

      if (c1 == 'T' .or. c1 == 't') then
         db_model(i)%has_db_energy = .true.
      else if (c1 == 'F' .or. c1 == 'f') then
         db_model(i)%has_db_energy = .false.
      else
         if (rangml == 0) write (6, *) 'ML: Error in reading T/F energy in database inp file'
         call mld_mpi_finalize('error 1 in reading db_file.inp, subroutine read_db_file')
         stop
      end if

      if (c2 == 'T' .or. c2 == 't') then
         db_model(i)%has_db_force = .true.
      else if (c2 == 'F' .or. c2 == 'f') then
         db_model(i)%has_db_force = .false.
      else
         if (rangml == 0) write (6, *) 'ML: Error in reading T/F force in database inp file'
         call mld_mpi_finalize('error 2 in reading db_file.inp, subroutine read_db_file')
         stop
      end if

      if (c3 == 'T' .or. c3 == 't') then
         db_model(i)%has_db_stress = .true.
      else if (c3 == 'F' .or. c3 == 'f') then
         db_model(i)%has_db_stress = .false.
      else
         if (rangml == 0) write (6, *) 'ML: Error in reading T/F  stress in database inp file'
         call mld_mpi_finalize('error 3 in reading db_file.inp, subroutine read_db_file')
         stop
      end if

      !debug write (*,*) db_model(i)%class, db_model(i)%klm, db_model(i)%no_total, db_model(i)%no_selec, db_model(i)%has_db_energy, db_model(i)%has_db_force, db_model(i)%has_db_stress, db_model(i)%w_e, db_model(i)%w_f, db_model(i)%w_s
   END DO

   rewind (inp)

   if (i .ne. n_db_class) then
      if (rangml == 0) write (6, *) 'Problems in readind db_file', inp
      if (rangml == 0) write (6, *) 'These values should be equal', n_db_class, i
      call mld_mpi_finalize("read_db_file")
      stop
   end if

   if (debug) then
      if (rangml == 0) write (6, *) 'ML: Reading DB file class, klm , no_total, no_selec'
      do i = 1, n_db_class
         if (rangml == 0) write (*, '(i3, "  ",(a)," ",(a),i6,i6)') i, db_model(i)%class, db_model(i)%klm, db_model(i)%no_total, db_model(i)%no_selec
         if (db_model(i)%no_total .lt. db_model(i)%no_selec) then
            if (rangml == 0) then
               write (6, *) 'Some problems in readinf db_file', db_file
               write (6, *) 'In the class', db_model(i)%class, 'and the KLM', db_model(i)%klm
               write (6, *) 'The number of total files', db_model(i)%no_total, 'is lower  than the no of selected files ', db_model(i)%no_selec
            end if
            call mld_mpi_finalize("no of files in read_db_file")
            stop
         end if
      end do
   end if

   iconf_temp = 0
   do i = 1, size(db_model)
      iconf_temp = iconf_temp + db_model(i)%no_selec
   end do
   iconf_data_train = iconf_temp

   iconf_temp = 0
   do i = 1, size(db_model)
      iconf_temp = iconf_temp + db_model(i)%no_total - db_model(i)%no_start + 1
   end do
   iconf_data = iconf_temp
   iconf_data_test = iconf_data - iconf_data_train

   if (allocated(config_real)) deallocate (config_real); allocate (config_real(iconf_data))
   if (allocated(config_desc)) deallocate (config_desc); allocate (config_desc(iconf_data))
   if (allocated(db_train)) deallocate (db_train); allocate (db_train(iconf_data_train))
   if (allocated(db_test)) deallocate (db_test); allocate (db_test(iconf_data_test))

   _MLD_END_
end subroutine read_db_file



subroutine prepare_name_file_from_db
   ! Read some information in the database_file and set some vectors
   ! Input:
   !       selection_type
   !       db_model object
   ! Output:
   !       db_train(1:iconf_data_train), db_test(1:iconf_data_test) - map the corresponding configuration between 1 to iconf_data
   !       config_real%
   !                  %class %klm %cnumber
   !                  %filename %i_order %train

   use ml_in_ndm_module, only: rangml, debug, &
      seed, train_only
   use module_db_setup, only: selection_type, selection_type_first, selection_type_last, &
      selection_type_random, selection_type_first_start, &
      db_path,  db_train, db_test, iconf_data, iconf_data_train, iconf_data_test
   use module_optimization, only: optimize_weights_db
   use derived_types, only: db_model, config_real
   use module_optimization, only: itopt, optimize_weights_db, optimize_weights_chem
   use module_covariance, only: train_covariance_matrix
   use mld_logger

   implicit none

   integer  :: i, i_conf, number, fnumber
   character(len=180)   :: ftemp, hdf5_ftemp, hdf5_gtemp
   character(len=7)     :: cnumber
   integer, dimension(:), allocatable     :: anum
   logical, dimension(:), allocatable     :: ltmp_sel
   integer  :: icount_temp, ii, i_train, i_test
   ! number of config to be analyzed
   logical  :: is_already_selected, lprint_local
   integer :: unumber1, unumber2

   _NAMECURRENT_("prepare_name_file_from_db")



   _MLD_BEGIN_
   number = 1e+6
   icount_temp = 0
   if (itopt == 0 ) then
      lprint_local = .true.
      if (train_covariance_matrix) lprint_local = .false.
      if (optimize_weights_chem) lprint_local = .false.
      if (optimize_weights_db) lprint_local = .false.
   else
      lprint_local = .false.
   end if

   if (lprint_local) then
      if (rangml == 0) write (*, '("ML:|---------------------------- Review of the DB model -----------------------------|")')
      if (rangml == 0) write (*, '("ML:  i_db_line    class klm           start          selected            total ")')
   end if
   ! preparing the name of the selected files for the traininf files
   config_real(:)%train = .false.
   do i = 1, size(db_model)
      if (lprint_local) then
         if (rangml == 0) write (6, '("ML:  ", i4, "           ", (a),"  ", (a), "        ", i6, "           ", i6,"              ", i6)') &
            i, db_model(i)%class, db_model(i)%klm, db_model(i)%no_start, db_model(i)%no_selec, db_model(i)%no_total
      end if
      if (allocated(ltmp_sel)) deallocate (ltmp_sel); allocate (ltmp_sel(db_model(i)%no_total))
      ltmp_sel(:) = .false.
      if (selection_type == selection_type_random) then
         if (allocated(anum)) deallocate (anum); allocate (anum(db_model(i)%no_selec))
         seed = seed + 1
         call rks2(db_model(i)%no_total, db_model(i)%no_selec, seed, anum)
      end if
      !i_start_name_file=10
      do i_conf = db_model(i)%no_start, db_model(i)%no_selec + db_model(i)%no_start - 1
         icount_temp = icount_temp + 1

         select case (selection_type)
          case (selection_type_first)
            fnumber = number + i_conf
          case (selection_type_last)
            fnumber = number + db_model(i)%no_total - db_model(i)%no_selec + i_conf
          case (selection_type_random)
            fnumber = number + anum(i_conf)
          case (selection_type_first_start)
            fnumber = number + i_conf
          case default
            if (rangml == 0) write (*, *) 'ML error: Only 4 possible choices for selection_type 1, 2 , 3 or 4. Read the manual'
            stop "fatal read_db_file"
         end select

         ! When drop_short_dist filtering compacted this class, fnumber-number
         ! is a position in the compacted (post-filter) sequence: remap it to
         ! the original on-disk POSCAR file number.
         if (allocated(db_model(i)%kept_file_number)) then
            fnumber = number + db_model(i)%kept_file_number(fnumber - number)
         end if

         ltmp_sel(i_conf) = .true.
         write (cnumber, '(i7)') fnumber
         ftemp = trim(adjustl(db_path))//db_model(i)%class//"_"//db_model(i)%klm//"_"//cnumber(2:7)//'.poscar'
         hdf5_ftemp = trim(adjustl(db_path))//db_model(i)%class//".h5"
         hdf5_gtemp = db_model(i)%class//"_"//db_model(i)%klm//"_"//cnumber(2:7)

         config_real(icount_temp)%class = db_model(i)%class
         config_real(icount_temp)%klm = db_model(i)%klm
         config_real(icount_temp)%cnumber = cnumber(2:7)
         config_real(icount_temp)%db_line = i
         config_real(icount_temp)%no_file_in_db_line = fnumber - number
         config_real(icount_temp)%filename = ftemp(1:80)
         config_real(icount_temp)%hdf5_filename = hdf5_ftemp(1:80)
         config_real(icount_temp)%hdf5_grpname = hdf5_gtemp(1:80)
         config_real(icount_temp)%train = .true.
         config_real(icount_temp)%selected = .true.
         !if (rangml==0) write (*,*) 'train', config_real(icount_temp)%filename, ltmp_sel(i_conf)
         config_real(icount_temp)%w_e = db_model(i)%w_e
         config_real(icount_temp)%w_f = db_model(i)%w_f
         config_real(icount_temp)%w_s = db_model(i)%w_s

         if (optimize_weights_db) then
            config_real(icount_temp)%w_e_end = db_model(i)%w_e_end
            config_real(icount_temp)%w_f_end = db_model(i)%w_f_end
            config_real(icount_temp)%w_s_end = db_model(i)%w_s_end
         end if

         config_real(icount_temp)%has_energy = db_model(i)%has_db_energy
         config_real(icount_temp)%has_force = db_model(i)%has_db_force
         config_real(icount_temp)%has_stress = db_model(i)%has_db_stress

         if (debug) then
            if (rangml == 0) write (*, '("ML: files of the database...in prepare_name_file_from_db: ",(a))') trim(ftemp)
         end if
         ! call read_poscar_sasha(ftemp,icount_temp)
      end do                  ! iconf

      !if (rangml==0) write (*,*) 'int 1', i, icount_temp
      if (db_model(i)%no_selec == db_model(i)%no_total) cycle

      do i_conf = db_model(i)%no_start, db_model(i)%no_total

         if (selection_type == selection_type_random) then
            is_already_selected = .false.
            do ii = 1, db_model(i)%no_selec
               if (i_conf == anum(ii)) is_already_selected = .true.
            end do
            if (is_already_selected) cycle
         else
            if (ltmp_sel(i_conf)) cycle
         end if
         !if (rangml==0) write (*,*) 'rest', config_real(icount_temp)%filename, config_real(icount_temp)%train
         icount_temp = icount_temp + 1
         fnumber = number + i_conf
         if (allocated(db_model(i)%kept_file_number)) then
            fnumber = number + db_model(i)%kept_file_number(fnumber - number)
         end if
         write (cnumber, '(i7)') fnumber
         ftemp = trim(adjustl(db_path))//db_model(i)%class//"_"//db_model(i)%klm//"_"//cnumber(2:7)//'.poscar'

         config_real(icount_temp)%class = db_model(i)%class
         config_real(icount_temp)%klm = db_model(i)%klm
         config_real(icount_temp)%cnumber = cnumber(2:7)
         config_real(icount_temp)%filename = ftemp(1:80)
         config_real(icount_temp)%train = .false.
         config_real(icount_temp)%db_line = i
         config_real(icount_temp)%no_file_in_db_line = fnumber - number
         config_real(icount_temp)%w_e = db_model(i)%w_e
         config_real(icount_temp)%w_f = db_model(i)%w_f
         config_real(icount_temp)%w_s = db_model(i)%w_s

         if (optimize_weights_db) then
            config_real(icount_temp)%w_e_end = db_model(i)%w_e_end
            config_real(icount_temp)%w_f_end = db_model(i)%w_f_end
            config_real(icount_temp)%w_s_end = db_model(i)%w_s_end
         end if

         config_real(icount_temp)%has_energy = db_model(i)%has_db_energy
         config_real(icount_temp)%has_force = db_model(i)%has_db_force
         config_real(icount_temp)%has_stress = db_model(i)%has_db_stress

      end do                  ! i_conf

      !if (rangml==0) write (*,*) 'int 2', i, icount_temp
   end do                  ! i size(db_model)

   if (icount_temp /= iconf_data) then
      if (rangml == 0) write (6, '("ML: Problems in the reading database icount_temp, iconf_data: ", 2i9)') &
         icount_temp, iconf_data
      stop "inconsistencies in the DB reading: subroutine prepare_name_file_from_db"
   end if

   if (lprint_local) then
      if (rangml == 0) write (*, '("ML:|---------------------------|")')
   end if

   !Fill the db_test and db_train vectors.
   i_train = 0
   i_test = 0
   do i = 1, iconf_data
      if (config_real(i)%train) then
         i_train = i_train + 1
         db_train(i_train) = i
      else
         i_test = i_test + 1
         db_test(i_test) = i
      end if
   end do

   if (i_train /= iconf_data_train) then
      if (rangml == 0) write (6, '("ML: Problems in the reading database icount_temp, iconf_data_train:  ")') i_train, iconf_data_train
      stop "inconsistencies in the DB train reading: subroutine prepare_name_file_from_db"
   end if
   if (i_test /= iconf_data_test) then
      if (rangml == 0) write (6, '("ML: Problems in the reading database icount_temp, iconf_data_test:  ")') i_test, iconf_data_test
      stop "inconsistencies in the DB test reading: subroutine prepare_name_file_from_db"
   end if

   if (rangml == 0) then
      open (file='train_files.milady', newunit=unumber1, action='write')
      if (.not. (train_only)) open (file='test_files.milady', newunit=unumber2, action='write')
      do i = 1, iconf_data
         if (config_real(i)%train) then
            write (unumber1, '(a80,i9)') trim(config_real(i)%filename), i
         else
            if ( .not. train_only)  write (unumber2, '(a80,i9)') trim(config_real(i)%filename), i
         end if
      end do
      close (unumber1, status='keep')
      if (.not. train_only)  close (unumber2, status='keep')
   end if

   _MLD_END_
end subroutine prepare_name_file_from_db


#if(MLD_HDF5)
subroutine fill_config_weights_and_mass_from_hdf5_file(ifile, nb_elements, elements)

   use derived_types, only: config_real
   use module_chemical_species, only: periodic_table_element, size_periodic_table, fix_no_of_elements, fix_ch_elements
   use ml_in_ndm_module, only:   fix_weighted_for_element, fix_weighted_for_element_3ch, &
      fix_weighted_for_element_ini, fix_weighted_for_element_3ch_ini, &
      weighted, weighted_3ch, &
      linvisible, fix_ch_elements_invisible, fix_no_of_elements_invisible
   use module_db_poscar, only: fix_ref_energy_per_element
   use mld_logger

   implicit none
   integer, intent(in)  :: ifile, nb_elements
   character(len=2), dimension(nb_elements), intent(in) :: elements
   integer :: i,j, icnt

   _NAMECURRENT_("fill_config_weights_and_mass_from_hdf5_file")

   _MLD_BEGIN_

   if (allocated(config_real(ifile)%ref_energy_per_element)) deallocate (config_real(ifile)%ref_energy_per_element); allocate (config_real(ifile)%ref_energy_per_element(nb_elements))
   if (allocated(config_real(ifile)%mass_per_type)) deallocate (config_real(ifile)%mass_per_type); allocate (config_real(ifile)%mass_per_type(nb_elements))
   if (allocated(config_real(ifile)%Z_per_type)) deallocate (config_real(ifile)%Z_per_type); allocate (config_real(ifile)%Z_per_type(nb_elements))
   if (allocated(config_real(ifile)%covalent_radius_per_type)) deallocate (config_real(ifile)%covalent_radius_per_type); allocate (config_real(ifile)%covalent_radius_per_type(nb_elements))
   if (allocated(config_real(ifile)%fix_type_poscar_to_periodic)) deallocate (config_real(ifile)%fix_type_poscar_to_periodic); allocate (config_real(ifile)%fix_type_poscar_to_periodic(nb_elements))
   if (allocated(config_real(ifile)%itype_to_global)) deallocate (config_real(ifile)%itype_to_global); allocate (config_real(ifile)%itype_to_global(nb_elements))

   ! assign the appropiate mass
   do i = 1, nb_elements
      do j = 1, size_periodic_table
         if (elements(i) == periodic_table_element(j)%symbol) then
            config_real(ifile)%fix_type_poscar_to_periodic(i) = periodic_table_element(j)%Z
            config_real(ifile)%mass_per_type(i) = periodic_table_element(j)%mass
            config_real(ifile)%Z_per_type(i) = periodic_table_element(j)%Z
            config_real(ifile)%covalent_radius_per_type(i) = periodic_table_element(j)%covalent_radius
         end if
      end do
   end do

   ! assign the appropiate weigths
   if (weighted) then
      do i = 1, nb_elements
         icnt = 0
         do j = 1, fix_no_of_elements
            if (fix_ch_elements(j) == elements(i)) then
               config_real(ifile)%itype_to_global(i) = j
               !$! if (weighted_auto) then
               !$!    fix_weighted_for_element(j) = dsqrt(config_real(ifile)%mass_per_type(i))
               !$!    if (weighted_3ch) then
               !$!       fix_weighted_for_element_3ch(j) = config_real(ifile)%Z_per_type(i)
               !$!    end if
               !$! end if
               icnt = icnt + 1
               config_real(ifile)%ref_energy_per_element(i) = fix_ref_energy_per_element(j)
            end if
         end do
         if (icnt == 0) then
            write (6, *) "ML: We have detected in the database an element that is not in the input list:"
            write (6, *) "ML: The unknown element is ", elements(i)
            write (6, *) "ML: Please update the fix_no_of_elements and chemical_elements."
            stop 'fatal in read_poscar, unknown element'
         end if
      end do
      if (linvisible) then
         if (allocated(config_real(ifile)%invisible_per_type)) deallocate (config_real(ifile)%invisible_per_type); allocate (config_real(ifile)%invisible_per_type(nb_elements))
         config_real(ifile)%invisible_per_type(:) = .false.
         do i = 1, nb_elements
            icnt = 0
            do j = 1, fix_no_of_elements_invisible
               if (fix_ch_elements_invisible(j) == elements(i)) then
                  config_real(ifile)%invisible_per_type(i) = .true.
               end if
            end do
         end do
      end if
   else                    ! not weighted
      do i = 1, nb_elements
         icnt = 0
         do j = 1, fix_no_of_elements
            if (fix_ch_elements(j) == elements(i)) then
               icnt = icnt + 1
               config_real(ifile)%ref_energy_per_element(i) = fix_ref_energy_per_element(j)
               config_real(ifile)%itype_to_global(i) = j
            end if
         end do
      end do
      if (icnt == 0) then
         write (6, *) "ML: 2 We have detected in the database an element that is not in the input list:"
         write (6, *) "ML: 2 The unknown element is ", elements(i)
         write (6, *) "ML: 2 Please update the fix_no_of_elements and chemical_elements."
         stop '2 fatal in read_poscar_sasha, unknown element'
      end if

   end if                  ! weighted
   !TODO maybe useless 
   if (weighted) then
      fix_weighted_for_element_ini(:) = fix_weighted_for_element(:)
      if (weighted_3ch) then
         fix_weighted_for_element_3ch_ini(:) = fix_weighted_for_element_3ch(:)
      end if
   end if

   _MLD_END_

end subroutine
#endif 


#if(MLD_HDF5)
subroutine fill_config_from_hdf5_file(ifile, E_fit1, E_fit2, E_total, has_energy, has_force, has_stress, nb_elements, elements, nspecies, im_local, box, xc, fp, st)

   use derived_types, only: config_real
   use module_chemical_species, only: size_periodic_table
   use ml_in_ndm_module, only: rangml, debug, im, imm
   !use module_db_poscar, only: fix_ref_energy_per_element
   use math
   use mld_logger

   implicit none
   integer, intent(in)  :: ifile, nb_elements, im_local
   integer, dimension(nb_elements), intent(in)    :: nspecies
   character(len=2), dimension(nb_elements), intent(in) :: elements
   real(kind=kind(0.d0)), dimension(3, 3), intent(in)  :: box
   real(kind=kind(1.d0)), dimension(3, im_local), intent(inout)  :: xc, fp
   real(kind=kind(1.d0)), intent(in)     :: E_total, E_fit1, E_fit2
   logical, intent(in)  :: has_energy, has_force, has_stress
   real(kind=kind(0.d0)), dimension(6), intent(in) :: st

   integer, dimension(:), allocatable     ::  ityp
   integer :: i,j, icnt
   real(kind=kind(1.d0)), dimension(:, :), allocatable      :: xp
   real(kind=kind(0.d0)), dimension(3, 3) :: box_inv
   real(kind=kind(1.d0)) :: tmp_ene, volume

   _NAMECURRENT_("fill_config_from_hdf5_file")

   _MLD_BEGIN_

   !write(*,*) "fill_from_hdf5_file NUMBER OF ELEMENTS" , nb_elements
   !write(*,*) "fill_from_hdf5_file NBSPECIES" , nspecies

   call fill_config_weights_and_mass_from_hdf5_file(ifile, nb_elements, elements)

   ! Number of atoms in the simulation box
   if (im_local .GT. imm) then
      imm = im_local
      if (debug) then
         write (6, '("ML: WARNING the imm values was changed to a upper value of im in fill_config_from_hdf5_file ", i7)') imm
      end if
   end if
   im = im_local

   icnt = 1
   if (allocated(ityp)) deallocate (ityp); allocate (ityp(1:im_local))
   do i = 1, nb_elements
      do j = 1, nspecies(i)
         ityp(icnt) = i
         icnt = icnt + 1
      end do
   end do

   ! Atom real coordinates
   ! cartesian coordinates are read
   if (allocated(xp)) deallocate (xp); allocate (xp(1:3, 1:im_local))
   call matinv_gen(box, box_inv)
   xp(1:3, 1:im_local) = xc(1:3, 1:im_local)
   xc(1:3, 1:im_local) = MatMul(box_inv(1:3, 1:3), xp(1:3, 1:im_local))

   if (allocated(config_real(ifile)%itype)) deallocate (config_real(ifile)%itype); allocate (config_real(ifile)%itype(im_local))
   if (allocated(config_real(ifile)%itype_db)) deallocate (config_real(ifile)%itype_db); allocate (config_real(ifile)%itype_db(im_local))
   if (allocated(config_real(ifile)%pos_cart)) deallocate (config_real(ifile)%pos_cart); allocate (config_real(ifile)%pos_cart(3, im_local))
   if (allocated(config_real(ifile)%pos_crst)) deallocate (config_real(ifile)%pos_crst); allocate (config_real(ifile)%pos_crst(3, im_local))
   if (allocated(config_real(ifile)%force)) deallocate (config_real(ifile)%force); allocate (config_real(ifile)%force(3, im_local))
   if (allocated(config_real(ifile)%atomic_spin)) deallocate (config_real(ifile)%atomic_spin); allocate (config_real(ifile)%atomic_spin(3, im_local))

   config_real(ifile)%ntypes = nb_elements
   config_real(ifile)%nat = im_local
   config_real(ifile)%itype(:) = ityp
   do i = 1, im_local
      config_real(ifile)%itype_db(i) = config_real(ifile)%itype_to_global(ityp(i))
   end do

   config_real(ifile)%pos_crst(1:3, 1:im_local) = xc(1:3, 1:im_local)
   config_real(ifile)%pos_cart(1:3, 1:im_local) = xp(1:3, 1:im_local)
   config_real(ifile)%force = fp

   tmp_ene = 0.d0
   do i = 1, im_local
      tmp_ene = tmp_ene + config_real(ifile)%ref_energy_per_element(ityp(i))
   end do

   config_real(ifile)%ref_energy = tmp_ene
   config_real(ifile)%energy(1) = E_total - tmp_ene
   config_real(ifile)%energy(2) = E_fit1 - tmp_ene
   config_real(ifile)%energy(3) = E_fit2 - tmp_ene

   config_real(ifile)%cell = box
   config_real(ifile)%stress(1:6) = st(1:6)
   config_real(ifile)%spin = 0

   config_real(ifile)%has_energy = has_energy .and. config_real(ifile)%has_energy
   config_real(ifile)%has_force = has_force .and. config_real(ifile)%has_force
   config_real(ifile)%has_stress = has_stress .and. config_real(ifile)%has_stress

   call calc_volume(config_real(ifile)%cell(:, 1), config_real(ifile)%cell(:, 2), config_real(ifile)%cell(:, 3), volume)
   config_real(ifile)%volume = volume

   ! write(*,*) "HDF FILLING"
   ! write(*,*) config_real(ifile)%ntypes, config_real(ifile)%nat
   ! write(*,*) config_real(ifile)%itype(:)
   ! write(*,*) config_real(ifile)%itype_db
   ! write(*,*) "!!!!!!!!!!!!!!!!!"
   ! write(*,*) config_real(ifile)%pos_crst
   ! write(*,*) "!!!!!!!!!!!!!!!!!"
   ! write(*,*) box_inv
   ! write(*,*) config_real(ifile)%pos_cart
   ! write(*,*) config_real(ifile)%force
   ! write(*,*) config_real(ifile)%ref_energy, config_real(ifile)%energy
   ! write(*,*) config_real(ifile)%cell
   ! write(*,*) config_real(ifile)%stress
   ! write(*,*) config_real(ifile)%spin
   ! write(*,*) config_real(ifile)%has_energy, config_real(ifile)%has_force, config_real(ifile)%has_stress
   ! write(*,*) config_real(ifile)%volume

   deallocate(xp); deallocate(ityp)

   _MLD_END_

end subroutine
#endif 

#if(MLD_HDF5)
subroutine read_hdf5_file(ifile)
   ! read hdf5 file containing multiple poscar files
   ! Input:
   !          ifile
   ! Output:
   !         im, imm
   !         imm is resized if im > imm
   !         config_real(ifile)%
   !                           %pos_cart  %pos_crst
   !                           %cell, %nat, %itype, %ntypes
   !                           %force %stress %spin
   !                           has_force, has_energy, has_stress
   use module_db_poscar, only: i_start_conf, i_final_conf
   use derived_types, only: config_real
   ! use module_db_setup, only: db_path
   use mld_logger
   use mld_hdf5, only: hFile, hdf5_open_file, hdf5_close_file, hdf5_open_group, &
                       hdf5_get_filename, hdf5_get_integer_group_attribute, &
                       hdf5_get_real_group_attribute, hdf5_get_string_group_dataset, &
                       hdf5_get_integer_group_dataset, hdf5_get_real_group_2Ddataset, &
                       hdf5_get_real_group_dataset

   implicit none

   integer, intent(in) :: ifile
   character(len=80)    :: name_file, name_group, hdf_fname
   integer  :: int_has_energy, int_has_force, int_has_stress
   logical :: has_energy, has_force, has_stress
   integer ::  nb_elements, im_local
   integer, dimension(:), allocatable     :: elements_mass, nspecies
   character(len=2), dimension(:), allocatable :: elements
   real(kind=kind(1.d0))      :: E_total, E_fit1, E_fit2, alat
   real(kind=kind(0.d0)), dimension(3, 3) :: box
   real(kind=kind(1.d0)), dimension(:, :), allocatable :: xc, fp
   real(kind=kind(1.d0)), dimension(6) :: st
   integer(8), parameter :: itwo =2 

   _NAMECURRENT_("read_hdf5_file")

   _MLD_BEGIN_

   ! opening file only if it's the first file we have to treat
   ! or if we need to read a new one (in which case, we also have to close the previous file we were reading)
   name_file = config_real(ifile)%hdf5_filename
   if(ifile == i_start_conf) then
      call hdf5_open_file(hFile,name_file)
   else
      call hdf5_get_filename(hFile, hdf_fname)
      if (name_file /= hdf_fname) then
         call hdf5_close_file(hFile)
         call hdf5_open_file(hFile,name_file)
      end if
   end if

   ! opening group
   name_group = config_real(ifile)%hdf5_grpname
   call hdf5_open_group(hFile,name_group)

   ! getting attributes
   call hdf5_get_integer_group_attribute(hFile,'nb_elements', nb_elements)
   if (nb_elements < 1) then
      print *, 'Number of chemical species can not be less than 1. Check yours DB files'
      stop "fatal hdf5 read"
   end if

   call hdf5_get_integer_group_attribute(hFile,'has_energy', int_has_energy)
   if (.not. (int_has_energy == 1 ) .or. (int_has_energy == 0 ) ) then 
     call log_warning(NAMECURRENT//" in hdf5 int_has_energy has illegal value "//vtoa(int_has_energy) )
   end if 

   if (int_has_energy == 1) then
      has_energy = .true.
   else
      has_energy = .false.
   end if
   call hdf5_get_integer_group_attribute(hFile,'has_force', int_has_force)

   if (.not. (int_has_force == 1 ) .or. (int_has_force == 0 ) ) then 
     call log_warning(NAMECURRENT//" in hdf5 int_has_force has illegal value "//vtoa(int_has_force) )
   end if 

   if (int_has_force == 1) then
      has_force = .true.
   else
      has_force = .false.
   end if
   call hdf5_get_integer_group_attribute(hFile,'has_stress', int_has_stress)

   if (.not. (int_has_stress == 1 ) .or. (int_has_stress == 0 ) ) then 
     call log_warning(NAMECURRENT//" in hdf5 int_has_stress has illegal value "//vtoa(int_has_stress) )
   end if 


   if (int_has_stress == 1) then
      has_stress = .true.
   else
      has_stress = .false.
   end if
   call hdf5_get_real_group_attribute(hFile,'Etot',E_total)
   call hdf5_get_real_group_attribute(hFile,'Efit1',E_fit1)
   call hdf5_get_real_group_attribute(hFile,'Efit2',E_fit2)
   call hdf5_get_real_group_attribute(hFile,'lat',alat)

   ! getting datasets
   if (allocated(elements)) deallocate (elements); allocate (elements(nb_elements))
   call hdf5_get_string_group_dataset(hFile, 'elements', nb_elements, itwo , elements)

   ! ToDo: mass is useless?????
   if (allocated(elements_mass)) deallocate (elements_mass); allocate (elements_mass(nb_elements))
   call hdf5_get_integer_group_dataset(hFile, 'mass', nb_elements, elements_mass)

   if (allocated(nspecies)) deallocate (nspecies); allocate (nspecies(nb_elements))
   call hdf5_get_integer_group_dataset(hFile, 'nb_species', nb_elements, nspecies)
   im_local = SUM(nspecies(:))

   call hdf5_get_real_group_2Ddataset(hFile, 'box', 3, 3, box)
   box(:, :) = alat*box(:, :)

   if (allocated(xc)) deallocate (xc); allocate (xc(1:3, 1:im_local))
   call hdf5_get_real_group_2Ddataset(hFile, 'xc', 3, im_local, xc)

   if (allocated(fp)) deallocate (fp); allocate (fp(1:3, 1:im_local))
   call hdf5_get_real_group_2Ddataset(hFile, 'fp', 3, im_local, fp)

   call hdf5_get_real_group_dataset(hFile, 'st', 6, st)

   call fill_config_from_hdf5_file(ifile, E_fit1, E_fit2, E_total, has_energy, has_force, has_stress, nb_elements, elements, nspecies, im_local, box, xc, fp, st)

   deallocate (elements); deallocate (elements_mass); deallocate(nspecies)
   deallocate(xc); deallocate(fp)

   if(ifile == i_final_conf) call hdf5_close_file(hFile)

   _MLD_END_

end subroutine
#endif 

subroutine read_poscar_sasha(ifile)
   ! read the Milady's poscar defined by Sasha
   ! Input:
   !          ifile, name_file
   ! Output:
   !         im, imm
   !         imm is resized if im > imm
   !         config_real(ifile)%
   !                           %pos_cart  %pos_crst
   !                           %cell, %nat, %itype, %ntypes
   !                           %force %stress %spin
   !                           has_force, has_energy, has_stress

   use ml_in_ndm_module, only: rangml, debug, im, imm, descriptor_type, &
      desc_forces, &
      fix_weighted_for_element, fix_weighted_for_element_3ch, &
      fix_weighted_for_element_ini, fix_weighted_for_element_3ch_ini, &
      weighted, weighted_3ch, &
      linvisible, fix_ch_elements_invisible, fix_no_of_elements_invisible
   use module_chemical_species, only: periodic_table_element, size_periodic_table, fix_ch_elements, &
      fix_no_of_elements
   use derived_types, only: config_real
   use module_db_poscar, only: fix_ref_energy_per_element
   use math
   use my_mpi_subroutines, only: my_broadcast_char, my_broadcast_val_int, my_broadcast_val_real, &
      my_broadcast_vect_real, my_broadcast_vect_int
   use mld_unit
   use mld_logger
   use mld_string
   use mld_mpi
   use mld_subworld, only: subworld, subrank

   implicit none

   character(len=80)    :: name_file
   integer, intent(in)  :: ifile
   logical  :: ok
   integer  :: inp
   integer  :: nb_elements
   character(len=3)     :: EFS_tag
   character(len=2), dimension(:), allocatable  :: element_poscar
   character(len=2), dimension(:), allocatable  :: internal_element_poscar
   character(len=2) :: ch_buffer
   real(kind_double), dimension(:), allocatable  :: vect_buffer
   integer, dimension(:), allocatable     :: mass_poscar
   real(kind=kind(1.d0))      :: E_total, E_fit1, E_fit2
   logical  :: has_energy, has_force, has_stress
   real(kind=kind(1.d0))      :: alat
   real(kind=kind(0.d0)), dimension(3, 3) :: box, box_inv
   character(len=80)    :: dummy, trimdummy
   character(len=500)   :: line
   integer  :: im_local, icnt, ij, ii, i_p, spin, i_pos_form, i, nitype, nitems2, itest_there_is_a_number, ios
   integer, dimension(:), allocatable     :: nspecies, ityp
   real(kind=kind(1.d0)), dimension(:, :), allocatable      :: xp, xc, fp, l_spin
   real(kind=kind(1.d0))      :: st(6), volume, tmp_ene

   _NAMECURRENT_("read_poscar_sasha")

   _MLD_BEGIN_

   name_file = config_real(ifile)%filename

   !openning ....
   if(subrank == 0) then
      inquire (file=trim(adjustl(name_file)), exist=ok)
      if (ok) then
         open (file=trim(adjustl(name_file)), newunit=inp, action='read')
      else
         write (6, *) 'ML: Fatal Error the poscar file is not there. The file: ', trim(adjustl(name_file)), rangml
         stop "fatal read_poscar_sasha"
      end if

      if (debug) then
         write (6, '("ML: Reading ..... in read_poscar_sasha ",(i5,a))') rangml, name_file
         write (6, '("ML: Reading ..... in read_poscar_sasha ",(i5, i6))') rangml, ifile
      end if
   end if
   !end_openning

   has_energy = .false.
   has_force = .false.
   has_stress = .false.
   if(subrank == 0) read (inp, *) EFS_tag, nb_elements
   call my_broadcast_char(EFS_tag, 0, subworld, codeml)
   call my_broadcast_val_int(nb_elements, 0, subworld, codeml)

   if (EFS_tag(1:1) == '1') has_energy = .true.
   if (EFS_tag(2:2) == '1') has_force = .true.
   if (EFS_tag(3:3) == '1') has_stress = .true.

   if (debug .AND. subrank == 0) then
      print *, 'DEBUG: found EFS', rangml, EFS_tag
      print *, 'DEBUG:', rangml, has_energy, has_force, has_stress
      print *, 'DEBUG: found', rangml, nb_elements, 'chemical elements'
   end if
   if(subrank == 0) BACKSPACE (inp)

   if (allocated(element_poscar)) deallocate (element_poscar); allocate (element_poscar(nb_elements))
   if (allocated(mass_poscar)) deallocate (mass_poscar); allocate (mass_poscar(nb_elements))

   ! IFORT_STYLE read  (inp,*) EFS_tag, nb_elements, ((element_poscar(ij),mass_poscar(ij)),ij=1,nb_elements), E_total, E_fit1, E_fit2

   if(subrank == 0) then
      read (inp, *, iostat=ios) EFS_tag, nb_elements, (element_poscar(ij), mass_poscar(ij), ij=1, nb_elements), E_total, E_fit1, E_fit2
      if (ios /= 0) then
         write(6,*) 'ERROR: end-of-file during read of POSCAR header line'
         write(6,*) '  POSCAR file : ', trim(name_file)
         write(6,*) '  nb_elements = ', nb_elements, '  iostat = ', ios
         call mld_mpi_abort('FATAL: cannot parse header in POSCAR file')
      end if
   end if
   call my_broadcast_val_real(E_total, 0, subworld, codeml)
   call my_broadcast_val_real(E_fit1 , 0, subworld, codeml)
   call my_broadcast_val_real(E_fit2 , 0, subworld, codeml)
   call my_broadcast_vect_int(mass_poscar, 0, subworld, codeml)
   do ii = 1, nb_elements
      ch_buffer = element_poscar(ii)
      call my_broadcast_char(ch_buffer, 0, subworld, codeml)
      element_poscar(ii) = trim(ch_buffer)
   end do

   if (debug) then
      ! IFORT_STYLE   if (rangml==0) write (6,*) EFS_tag, nb_elements, ((element_poscar(ij),mass_poscar(ij)), ij=1,nb_elements), E_total, E_fit1, E_fit2
   end if
   ! allocate(config_real(ifile)%weight_per_type(nb_elements))
   if (allocated(config_real(ifile)%ref_energy_per_element)) deallocate (config_real(ifile)%ref_energy_per_element); allocate (config_real(ifile)%ref_energy_per_element(nb_elements))
   if (allocated(config_real(ifile)%mass_per_type)) deallocate (config_real(ifile)%mass_per_type); allocate (config_real(ifile)%mass_per_type(nb_elements))
   if (allocated(config_real(ifile)%Z_per_type)) deallocate (config_real(ifile)%Z_per_type); allocate (config_real(ifile)%Z_per_type(nb_elements))
   if (allocated(config_real(ifile)%covalent_radius_per_type)) deallocate (config_real(ifile)%covalent_radius_per_type); allocate (config_real(ifile)%covalent_radius_per_type(nb_elements))
   if (allocated(config_real(ifile)%fix_type_poscar_to_periodic)) deallocate (config_real(ifile)%fix_type_poscar_to_periodic); allocate (config_real(ifile)%fix_type_poscar_to_periodic(nb_elements))
   if (allocated(config_real(ifile)%itype_to_global)) deallocate (config_real(ifile)%itype_to_global); allocate (config_real(ifile)%itype_to_global(nb_elements))


   ! assign the appropiate mass
   do ij = 1, nb_elements
      do i_p = 1, size_periodic_table
         if (element_poscar(ij) == periodic_table_element(i_p)%symbol) then
            config_real(ifile)%fix_type_poscar_to_periodic(ij) = periodic_table_element(i_p)%Z
            config_real(ifile)%mass_per_type(ij) = periodic_table_element(i_p)%mass
            config_real(ifile)%Z_per_type(ij) = periodic_table_element(i_p)%Z
            config_real(ifile)%covalent_radius_per_type(ij) = periodic_table_element(i_p)%covalent_radius
         end if
      end do
   end do

   ! assign the appropiate weigths
   if (weighted) then

      do ij = 1, nb_elements
         icnt = 0
         do i_p = 1, fix_no_of_elements
            if (fix_ch_elements(i_p) == element_poscar(ij)) then
               config_real(ifile)%itype_to_global(ij) = i_p
               !$! if (weighted_auto) then
               !$!    fix_weighted_for_element(i_p) = dsqrt(config_real(ifile)%mass_per_type(ij))
               !$!    if (weighted_3ch) then
               !$!       fix_weighted_for_element_3ch(i_p) = config_real(ifile)%Z_per_type(ij)
               !$!    end if
               !$! else
               !$!    !dd! config_real(ifile)%weight_per_type(ij) = fix_weighted_for_element(i_p)
               !$!    !dd! if (weighted_3ch) config_real(ifile)%weight_per_type_3ch(ij) = fix_weighted_for_element_3ch(i_p)
               !$! end if
               icnt = icnt + 1
               config_real(ifile)%ref_energy_per_element(ij) = fix_ref_energy_per_element(i_p)
            end if
         end do
         if (icnt == 0) then
            write (6, *) "ML: We have detected in the database an element that is not in the input list:"
            write (6, *) "ML: The unknown element is ", element_poscar(ij), " in the file ", name_file
            write (6, *) "ML: Please update the fix_no_of_elements and chemical_elements."
            stop 'fatal in read_poscar_sasha, unknown element'
         end if
      end do
      if (linvisible) then
         if (allocated(config_real(ifile)%invisible_per_type)) deallocate (config_real(ifile)%invisible_per_type); allocate (config_real(ifile)%invisible_per_type(nb_elements))
         config_real(ifile)%invisible_per_type(:) = .false.
         do ij = 1, nb_elements
            icnt = 0
            do i_p = 1, fix_no_of_elements_invisible
               if (fix_ch_elements_invisible(i_p) == element_poscar(ij)) then
                  config_real(ifile)%invisible_per_type(ij) = .true.
               end if
            end do
         end do
      end if

   else                    ! not weighted

      do ij = 1, nb_elements
         icnt = 0
         do i_p = 1, fix_no_of_elements
            if (fix_ch_elements(i_p) == element_poscar(ij)) then
               icnt = icnt + 1
               config_real(ifile)%ref_energy_per_element(ij) = fix_ref_energy_per_element(i_p)
               config_real(ifile)%itype_to_global(ij) = i_p
            end if
         end do
      end do
      if (icnt == 0) then
         write (6, *) "ML: 2 We have detected in the database an element that is not in the input list:"
         write (6, *) "ML: 2 The unknown element is ", element_poscar(ij), " in the file ", name_file
         write (6, *) "ML: 2 Please update the fix_no_of_elements and chemical_elements."
         stop '2 fatal in read_poscar_sasha, unknown element'
      end if
      !debug write(*, * ) 'id', rangml,   ifile, config_real(ifile)%itype_to_global(:)

   end if                  ! weighted

   !TODO the next seems to be useless 
   if (weighted) then
      fix_weighted_for_element_ini(:) = fix_weighted_for_element(:)
      if (weighted_3ch) then
         fix_weighted_for_element_3ch_ini(:) = fix_weighted_for_element_3ch(:)
      end if
   end if

   if (debug) then
      if (subrank == 0 ) then 
      do ij = 1, nb_elements
         write (6, *) 'in poscar type ', config_real(ifile)%filename
      end do
      end if 
   end if

   if (nb_elements < 1) then
      print *, 'Number of chemical spieces can not be less than 1. Check yours DB files'
      stop "fatal read_poscar_sasha"
   end if

   ! Lattice vector coordinates of the cell
   !read_inp ....
   if(subrank == 0) then
      read (inp, *) alat
      read (inp, *) box(1:3, 1)
      read (inp, *) box(1:3, 2)
      read (inp, *) box(1:3, 3)
   end if
   !broadcast_inp ....
   call my_broadcast_val_real(alat, 0, subworld, codeml)
   if (allocated(vect_buffer)) deallocate(vect_buffer) ; allocate(vect_buffer(size(box,1)))
   do ii = 1, size(box,2)
      vect_buffer = box(:,ii)
      call my_broadcast_vect_real(vect_buffer, 0, subworld, codeml)
      box(:,ii) = vect_buffer(:)
   end do

   box(:, :) = alat*box(:, :)

   !d! write(*,*) 'tot',rangml, alat
   !d! write(*,*) 'moto', box(:,2)

   ! read_inp
   if(subrank == 0) read (inp, '(a)') line
   call my_broadcast_char(line, 0, subworld, codeml)
   !d! write(*,*) rangml, 'toto--|', line

   if (itest_there_is_a_number(line) == 0) then
      ! if (debug) then
      !   if (rangml==0) write (*,*) "ML: POSCAR type  VASP5 file"
      ! end if
      nitype = nitems2(line)
      if (allocated(internal_element_poscar)) deallocate(internal_element_poscar)
      allocate(internal_element_poscar(nb_elements))
      read(line,*, iostat=ios) internal_element_poscar(:)
      if (ios /= 0) then
         write(6,*) 'ERROR: end-of-file during internal read of species line'
         write(6,*) '  POSCAR file : ', trim(name_file)
         write(6,*) '  line content: |', trim(line), '|'
         write(6,*) '  nb_elements = ', nb_elements, '  nitype = ', nitype, '  iostat = ', ios
         call mld_mpi_abort('FATAL: cannot parse species line in POSCAR header')
      end if
      do ii = 1, size(internal_element_poscar)
         if ( internal_element_poscar(ii) /= element_poscar(ii) ) then
            call log_warning('there is no consistency between the header of poscar file and VASP5 format') 
            call log_warning('  - the line with chemical symbol of species '//vtoa(ii))
            call log_warning('here are the species '//internal_element_poscar(ii)//'  '//element_poscar(ii))
            call mld_mpi_abort('FATAL ERROR: put the database in right format VASP5')
         end if
      end do
   else
      ! if (debug) then
      !   if (rangml==0) write (*,*) "ML: POSCAR type  VASP4 file"
      ! end if
      nitype = nitems2(line)
      if(subrank == 0) BACKSPACE (inp)
   end if

   if (nitype /= nb_elements) then
      write (6, *) 'There are inconsistencies in the number of species in the poscar file nb_elements, nitype', nb_elements, nitype
      write (6, *) 'Check this file:  ', name_file
      stop 'species wrong defined: subroutine read_poscar_sasha'
   end if

   if (allocated(nspecies)) deallocate (nspecies); allocate (nspecies(nitype))

   ! Number of atoms in the cell
   nspecies(:) = 0
   if(subrank == 0) then
      read (inp, *) (nspecies(i), i=1, nitype)
      read (inp, '(a)') dummy
   end if
   trimdummy = trim(adjustl(dummy))
   dummy=trimdummy
   call my_broadcast_vect_int(nspecies, 0, subworld, codeml)
   call my_broadcast_char(dummy, 0, subworld, codeml)

   !d! write(*,*) rangml, 'toto', nspecies
   !d! write(*,*) rangml, 'moto', dummy, len(dummy)

   trimdummy = trim(adjustl(dummy))
   if (.not. ((trimdummy(1:1) .eq. 'c') .or. (trimdummy(1:1) .eq. 'C') .or. (trimdummy(1:1) .eq. 'd') .or. (trimdummy(1:1) .eq. 'D'))) then
      write (6, *) 'The POSCAR file is not in the good format. The first letter in that line is not D , d, C or c'
      write (6, *) trimdummy(1:1), " ", dummy
      write (6, *) 'Check this file:  ', name_file
      stop "vasp format read_poscar_sasha"
   end if

   if ((trimdummy(1:1) .eq. 'c') .or. (trimdummy(1:1) .eq. 'C')) then
      i_pos_form = 0          ! cartesian format
   end if

   if ((trimdummy(1:1) .eq. 'd') .or. (trimdummy(1:1) .eq. 'D')) then
      i_pos_form = 1          ! crystalografic format
   end if

   ! Number of atoms in the simulation box
   im_local = SUM(nspecies(:))
   if (im_local .GT. imm) then
      imm = im_local
      if (debug) then
         write (6, '("ML: WARNING the imm values was changed to a upper value of im in read_poscar_sasha ", i7)') imm
      end if
   end if
   im = im_local

   ! write (*,*) im_local
   ! allocate(ityp(1:im_local))
   if (allocated(ityp)) deallocate (ityp); allocate (ityp(1:im_local))
   if (allocated(xc)) deallocate (xc); allocate (xc(1:3, 1:im_local))
   if (allocated(xp)) deallocate (xp); allocate (xp(1:3, 1:im_local))
   if (allocated(fp)) deallocate (fp); allocate (fp(1:3, 1:im_local))
   if (allocated(l_spin)) deallocate (l_spin); allocate (l_spin(1:3, 1:im_local))
   icnt = 0
   do ij = 1, nitype
      do i = 1, nspecies(ij)
         icnt = icnt + 1
         ityp(icnt) = ij
         ! read_inp
         if(subrank == 0) read (inp, *) xc(1:3, icnt)
      end do
   end do
   !broadcast_inp
   if (allocated(vect_buffer)) deallocate(vect_buffer) ; allocate(vect_buffer(size(xc,2)))
   do ii = 1, 3
      vect_buffer(:) = xc(ii,:)
      call my_broadcast_vect_real(vect_buffer, 0, subworld, codeml)
      xc(ii,:) = vect_buffer(:)
   end do

   ! Atom real coordinates
   if (i_pos_form == 1) then
      ! inthat case crystalografic coord are read
      xp(1:3, 1:im_local) = MatMul(box(1:3, 1:3), xc(1:3, 1:im_local))
   else if (i_pos_form == 0) then
      ! in that case cartesian coordinates are read
      call matinv_gen(box, box_inv)
      xp(1:3, 1:im_local) = xc(1:3, 1:im_local)
      xc(1:3, 1:im_local) = MatMul(box_inv(1:3, 1:3), xp(1:3, 1:im_local))
   end if
   if (desc_forces) then
      if (subrank == 0) then
         read (inp, *)
         do i = 1, im_local
            read (inp, *) fp(1:3, i)
         end do
      end if
      if (allocated(vect_buffer)) deallocate(vect_buffer) ; allocate(vect_buffer(size(fp,2)))
      do ii = 1,3
         vect_buffer(:) = fp(ii,:)
         call my_broadcast_vect_real(vect_buffer, 0, subworld, codeml)
         fp(ii,:) = vect_buffer(:)
      end do
      !d! write(*,*) rangml, 'toto', fp(:,10)

      if(subrank == 0) then
         read (inp, *)
         read (inp, *) st(1:6)
      end if
      call my_broadcast_vect_real(st, 0, subworld, codeml)
      !d! write(*,*) rangml, 'toto', st

      !stop "-- la vita e bella -- "

   end if

   !debug read  (inp,*)
   !debug read  (inp,*) spin
   spin = 0
   !call my_broadcast_val_int(spin, subrank, subworld, codeml) !AK: useless if we defined spin just above for everyone



   if (allocated(config_real(ifile)%itype)) deallocate (config_real(ifile)%itype); allocate (config_real(ifile)%itype(im_local))
   if (allocated(config_real(ifile)%itype_db)) deallocate (config_real(ifile)%itype_db); allocate (config_real(ifile)%itype_db(im_local))
   if (allocated(config_real(ifile)%pos_cart)) deallocate (config_real(ifile)%pos_cart); allocate (config_real(ifile)%pos_cart(3, im_local))
   if (allocated(config_real(ifile)%pos_crst)) deallocate (config_real(ifile)%pos_crst); allocate (config_real(ifile)%pos_crst(3, im_local))
   if (allocated(config_real(ifile)%force)) deallocate (config_real(ifile)%force); allocate (config_real(ifile)%force(3, im_local))
   if (allocated(config_real(ifile)%atomic_spin)) deallocate (config_real(ifile)%atomic_spin); allocate (config_real(ifile)%atomic_spin(3, im_local))

   config_real(ifile)%ntypes = nitype
   config_real(ifile)%nat = im_local
   config_real(ifile)%itype(:) = ityp

   do ii = 1, im_local
      config_real(ifile)%itype_db(ii) = config_real(ifile)%itype_to_global(ityp(ii))
   end do

   config_real(ifile)%pos_crst(1:3, 1:im_local) = xc(1:3, 1:im_local)
   config_real(ifile)%pos_cart(1:3, 1:im_local) = xp(1:3, 1:im_local)
   config_real(ifile)%force = fp

   tmp_ene = 0.d0
   do ii = 1, im_local
      tmp_ene = tmp_ene + config_real(ifile)%ref_energy_per_element(ityp(ii))
   end do

   config_real(ifile)%ref_energy = tmp_ene
   config_real(ifile)%energy(1) = E_total - tmp_ene
   config_real(ifile)%energy(2) = E_fit1 - tmp_ene
   config_real(ifile)%energy(3) = E_fit2 - tmp_ene

   config_real(ifile)%cell = box
   config_real(ifile)%stress(1:6) = st(1:6)
   config_real(ifile)%spin = spin
   config_real(ifile)%atomic_spin = l_spin


   config_real(ifile)%has_energy = has_energy .and. config_real(ifile)%has_energy
   config_real(ifile)%has_force = has_force .and. config_real(ifile)%has_force
   config_real(ifile)%has_stress = has_stress .and. config_real(ifile)%has_stress


   call calc_volume(config_real(ifile)%cell(:, 1), config_real(ifile)%cell(:, 2), config_real(ifile)%cell(:, 3), volume)
   config_real(ifile)%volume = volume

   if(subrank == 0) close (inp)

   ! if(ifile == 1) then
   !    write(*,*) "REF FILLING"
   !    write(*,*) config_real(ifile)%ntypes, config_real(ifile)%nat
   !    write(*,*) config_real(ifile)%itype(:)
   !    write(*,*) config_real(ifile)%itype_db
   !    write(*,*) "!!!!!!!!!!!!!!!!!"
   !    write(*,*) config_real(ifile)%pos_crst
   !    write(*,*) "!!!!!!!!!!!!!!!!!"
   !    write(*,*) box_inv
   !    write(*,*) config_real(ifile)%pos_cart
   !    write(*,*) config_real(ifile)%force
   !    write(*,*) config_real(ifile)%ref_energy, config_real(ifile)%energy
   !    write(*,*) config_real(ifile)%cell
   !    write(*,*) config_real(ifile)%stress
   !    write(*,*) config_real(ifile)%spin
   !    write(*,*) config_real(ifile)%has_energy, config_real(ifile)%has_force, config_real(ifile)%has_stress
   !    write(*,*) config_real(ifile)%volume
   ! end if

   _MLD_END_
end subroutine read_poscar_sasha


subroutine fix_poscar_weights
   use ml_in_ndm_module, only: rangml, weighted, weighted_auto, weighted_3ch, &
      renorm_mass, renorm_cov, &
      fix_weighted_for_element, fix_weighted_for_element_3ch, &
      fix_weighted_for_element_ini, fix_weighted_for_element_3ch_ini
   use module_chemical_species, only:    fix_mass_elements
   use mld_logger
   implicit none
   _NAMECURRENT_("fix_poscar_weights")

   _MLD_BEGIN_

   if (weighted) then
      if (weighted_auto) then
         if (size(fix_mass_elements) == 0) then
            if (rangml == 0) write (6, '("ML: fix_mass_elements has 0 componeents. Is not possible because weighted_auto and weighted is true")')
            stop 'fix_mass_elements error in prepare_train_dimensions'
         else
            !dnc renorm_mass = dsqrt(sum(fix_mass_elements(:))/size(fix_mass_elements(:)))
            renorm_mass = dsqrt(sum(fix_weighted_for_element_ini(:)**2)/size(fix_weighted_for_element_ini(:)))
         end if

         if ((weighted_3ch) .and. (weighted_auto)) then
            if (size(fix_weighted_for_element_3ch) == 0) then
               if (rangml == 0) write (6, '("ML: fix_Z elements has 0 components. Is not possible because weighted_auto and weighted_3ch is true")')
               stop 'fix_Z elements error in prepare_train_dimensions'
            else
               renorm_cov = sum(fix_weighted_for_element_3ch_ini(:))/size(fix_weighted_for_element_3ch_ini(:))
            end if
         end if
         fix_weighted_for_element(:) = fix_weighted_for_element_ini(:)/renorm_mass
         if (weighted_3ch) fix_weighted_for_element_3ch(:) = fix_weighted_for_element_3ch_ini(:)/renorm_cov
      else
         fix_weighted_for_element(:) = fix_weighted_for_element_ini(:)
         if (weighted_3ch) fix_weighted_for_element_3ch(:) = fix_weighted_for_element_3ch_ini(:)
      end if
   end if

   _MLD_END_

end subroutine fix_poscar_weights

subroutine fix_atoms_weights(iconf)
   use ml_in_ndm_module, only: weighted,  weighted_3ch, & 
      fix_weighted_for_element, fix_weighted_for_element_3ch
   use derived_types, only: config_real
   use module_chemical_species, only:    img_num_ch, img_weighted, fix_wspecies
   use mld_logger
   implicit none
   integer , intent(in) :: iconf
   integer :: ij, nb_elements, i_p, ich
   _NAMECURRENT_("fix_atoms_weights")

   _MLD_BEGIN_
   nb_elements = config_real(iconf)%ntypes
   if (allocated(config_real(iconf)%weight_per_type)) deallocate (config_real(iconf)%weight_per_type); allocate (config_real(iconf)%weight_per_type(nb_elements))
   if (allocated(config_real(iconf)%weight_per_type_3ch)) deallocate (config_real(iconf)%weight_per_type_3ch); allocate (config_real(iconf)%weight_per_type_3ch(nb_elements))

   if (weighted) then
      do ij = 1, nb_elements
         config_real(iconf)%weight_per_type(ij) = fix_weighted_for_element(config_real(iconf)%itype_to_global(ij))
         if (weighted_3ch) config_real(iconf)%weight_per_type_3ch(ij) =  fix_weighted_for_element_3ch(config_real(iconf)%itype_to_global(ij))
      end do
   else
      do ij = 1, nb_elements
         config_real(iconf)%weight_per_type(ij) = 1.d0
      end do
   end if


   ! image_arnaud
   if (img_weighted) then
      if (allocated(config_real(iconf)%wspecies_per_type_ch)) deallocate(config_real(iconf)%wspecies_per_type_ch)
      allocate(config_real(iconf)%wspecies_per_type_ch(nb_elements,img_num_ch))
   end if
   !end  image_arnaud

   ! image_arnaud
   if (img_weighted) then
      do ij = 1, nb_elements
         i_p =  config_real(iconf)%itype_to_global(ij)
         do ich = 1 , img_num_ch
            config_real(iconf)%wspecies_per_type_ch(ij,ich) = fix_wspecies(i_p,ich)
         end do
      end do
   end if
   ! end image_arnaud

   _MLD_END_

end subroutine fix_atoms_weights



subroutine combinations(n, k, n_sample, seedin, b)

   implicit none

   integer(kind=4), intent(in)      :: n, k, seedin, n_sample
   integer, dimension(n_sample, k), intent(out) :: b

   integer, dimension(k)      :: a
   integer(kind=4)      :: i, j, c, seed
   logical  :: test

   seed = seedin
   call rks2(n, k, seed, a)
   b(1, :) = a(:)
   i = 2
   c = 0
   do while (i - 1 .ne. n_sample)
      c = c + 1
      seed = seed + c
      call rks2(n, k, seed, a)
      test = .true.
      do j = 1, i - 1
         if (all(a(:) == b(j, :))) then
            test = .false.
            exit
         end if
      end do
      if (test) then
         b(i, :) = a(:)
         i = i + 1
      end if
   end do
end subroutine combinations


subroutine rks2(n, k, seed, a)

   use ml_in_ndm_module, only: rangml
   implicit none

   integer(kind=4)      :: k, c1, c2, i, k0, n, seed
   integer(kind=4), dimension(k)    :: a
   real(kind=8)   :: r, r8_uniform_01

   if (k < 0 .or. n < k) then
      if (rangml == 0) write (*, '(a)') ''
      if (rangml == 0) write (*, '(a)') 'KSUB_RANDOM2 - Fatal error!'
      if (rangml == 0) write (*, '(a,i8)') '  N = ', n
      if (rangml == 0) write (*, '(a,i8)') '  K = ', k
      if (rangml == 0) write (*, '(a)') '  but 0 <= K <= N is required!'
      stop 1
   end if

   if (k == 0) return
   c1 = k
   c2 = n
   k0 = 0
   i = 0
   do i = 1, n
      r = r8_uniform_01(seed)
      if (real(c2, kind=8)*r <= real(c1, kind=8)) then
         c1 = c1 - 1
         k0 = k0 + 1
         a(k0) = i
         if (c1 <= 0) then
            exit
         end if
      end if
      c2 = c2 - 1
   end do
end subroutine rks2



function r8_uniform_01(seed)
   ! R8_UNIFORM_01 returns a unit pseudorandom R8.
   !  Discussion:
   !    An R8 is a real ( kind = 8 ) value.
   !    For now, the input quantity SEED is an integer ( kind = 4 ) variable.
   !    This routine implements the recursion
   !      seed = 16807 * seed mod ( 2^1 - 1 )
   !      r8_uniform_01 = seed / ( 2^31 - 1 )
   !
   !    The integer arithmetic never requires more than 32 bits,
   !    including a sign bit.
   !
   !    If the initial seed is 12345, then the first three computations are
   !      Input     Output      R8_UNIFORM_01
   !      SEED      SEED
   !         12345   207482415  0.096616
   !     207482415  1790989824  0.833995
   !    1790989824  2035175616  0.947702
   !
   !  Licensing: This code is distributed under the GNU LGPL license.
   !  Modified: 05 July 2006
   !  Author: John Burkardt
   !
   !  Reference:
   !    Paul Bratley, Bennett Fox, Linus Schrage,
   !    A Guide to Simulation,
   !    Springer Verlag, pages 201-202, 1983.
   !
   !    Bennett Fox,
   !    Algorithm 647:
   !    Implementation and Relative Efficiency of Quasirandom
   !    Sequence Generators,
   !    ACM Transactions on Mathematical Software,
   !    Volume 12, Number 4, pages 362-376, 1986.
   !
   !    Pierre LEcuyer,
   !    Random Number Generation,
   !    in Handbook of Simulation,
   !    edited by Jerry Banks,
   !    Wiley Interscience, page 95, 1998.
   !
   !    Peter Lewis, Allen Goodman, James Miller
   !    A Pseudo-Random Number Generator for the System/360,
   !    IBM Systems Journal,
   !    Volume 8, pages 136-143, 1969.
   !
   !  Parameters:
   !    Input/output, integer ( kind = 4 ) SEED, the "seed" value, which should
   !    NOT be 0. On output, SEED has been updated.
   !    Output, real ( kind = 8 ) R8_UNIFORM_01, a new pseudorandom variate,
   !    strictly between 0 and 1.

   use ml_in_ndm_module, only: rangml
   implicit none

   integer(kind=4)      :: i4_huge, k, seed
   real(kind=8)   :: r8_uniform_01

   if (seed == 0) then
      if (rangml == 0) write (*, '(a)') ''
      if (rangml == 0) write (*, '(a)') 'R8_UNIFORM_01 - Fatal error!'
      if (rangml == 0) write (*, '(a)') '  Input value of SEED = 0.'
      stop 1
   end if

   k = seed/127773

   seed = 16807*(seed - k*127773) - k*2836

   if (seed < 0) then
      seed = seed + i4_huge()
   end if
   !
   !  Although SEED can be represented exactly as a 32 bit integer,
   !  it generally cannot be represented exactly as a 32 bit real number!
   !
   r8_uniform_01 = real(seed, kind=8)*4.656612875D-10
end function r8_uniform_01



function i4_huge()
   ! I4_HUGE returns a "huge" I4.
   !    This code is distributed under the GNU LGPL license.
   !    Output, integer ( kind = 4 ) I4_HUGE, a "huge" integer.
   implicit none
   integer(kind=4) i4_huge
   i4_huge = 2147483647
end function i4_huge


subroutine read_mask(ifile)
! this subroutine can be used only if lmask = . true.
#if (PARA)
   !use mpi
   use mld_mpi, only:  comm_mld
#endif
   use derived_types, only: config_real, config_desc
   use ml_in_ndm_module, only: rangml, mask_file
   use module_db_setup, only: db_path
   use mld_unit
   use mld_logger

   implicit none

   character(len=80)    :: name_file
   character(len=200) :: inputline
   integer, intent(in)  :: ifile
   logical  :: ok
   integer  :: unitmask, iiostat, numlines, ii,itmp


   _NAMECURRENT_("read_mask")

   _MLD_BEGIN_

   if (mask_file) then
      name_file = "mask_model.in"
   else
      name_file = trim(adjustl(db_path))//"/"//config_real(ifile)%class//"_"//config_real(ifile)%klm//"_"//config_real(ifile)%cnumber//".mask"
   end if
   name_file=trim(adjustl(name_file))

   inquire (file=trim(adjustl(name_file)), exist=ok)
   if (ok) then
      open (file=trim(adjustl(name_file)), newunit=unitmask, action='read')
   else
      if (rangml == 0) write (6, *) 'ML: Fatal Error the mask file is not there. The file: ', trim(adjustl(name_file))
      stop "fatal read_mask"
   end if

   if (allocated(config_desc(ifile)%amask)) deallocate(config_desc(ifile)%amask)
   allocate(config_desc(ifile)%amask(config_real(ifile)%nat))

   if (rangml ==0) then

      open(newunit=unitmask, file=name_file, status='old', action='read', position='rewind')
      numlines = 0
      loop1: do
         read(unitmask,*,iostat=iiostat) inputline
         if (iiostat < 0) then
            exit loop1
         end if
         numlines = numlines + 1
      end do loop1

      rewind(unitmask)

      config_desc(ifile)%amask(:) = .false.
      do ii = 1, numlines
         read(unitmask,*) itmp
         if (itmp > config_real(ifile)%nat) then
            call log_warning("saving the descriptor of this atom is not possible being out of box range")
         end if
         config_desc(ifile)%amask(itmp) = .true.
      end do

      close(unitmask)


   end if

   !TORC! call mpi_bcast(config_desc(ifile)%amask(:),config_real(ifile)%nat , MPI_LOGICAL, 0, MPI_COMM_WORLD, codeml)
   call comm_mld%bcast(0, config_desc(ifile)%amask(:))
   _MLD_END_

end subroutine read_mask
