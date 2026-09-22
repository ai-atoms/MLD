! HND XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
! HND X
! HND X   MiLaDy (Machine Learning Dynamics)
! HND X
! HND X   Milady was written and designed by Mihai-Cosmin Marinica, Alexandra M. Goryaeva
! HND X   Copyright 2015-2026.
! HND X
! HND X   MiLaDy is published and distributed under the
! HND X      Academic Software License v1.0 (ASL)
! HND X
! HND XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

#include "../MLD_MACROS.INC"

! ============================================================================
! MODULE module_json
!   Direct reader for large nested-JSON trajectory datasets, e.g. the
!   Materials Project MPtrj format:
!     { "mp-XXXXXXX": { "mp-XXXXXXX-N-M": {
!          "structure": {"lattice": {"matrix": [[..],[..],[..]]},
!                        "sites": [{"species":[{"element":"Fe"}], "xyz":[x,y,z]}, ...]},
!          "corrected_total_energy": ..., "force": [[fx,fy,fz],...] or null,
!          "stress": [[..],[..],[..]] or null, ... }, ... }, ... }
!   Activated when db_path ends with .json
!
!   No JSON library is used or linked in this project, and files of this
!   kind can be tens of GB with no newlines at all, so this module never
!   loads more than a small local window of the file into memory: a bulk
!   streaming pass finds every "structure" object's byte offset (marker
!   search over fixed-size overlapping chunks), then each configuration's
!   metadata/atoms are extracted on demand via small positioned reads and a
!   hand-rolled, string-aware brace/bracket scanner - this is not a general
!   JSON parser, only what is needed to extract MPtrj's fixed schema.
!
!   This module is entirely separate from the POSCAR and XYZ reading paths
!   (module_extxyz, read_poscar_sasha) - neither is touched by this file.
! ============================================================================
module module_json
   use module_kind_variables, only: kind_double
   implicit none
   private

   !> Whether the database is in MPtrj-style JSON format (auto-detected from db_path)
   logical, public :: db_json = .false.

   !> Bulk marker-search chunk size (bytes) and inter-chunk overlap (must
   !> exceed the marker length so a boundary-split marker is never missed)
   integer, parameter :: CHUNK_SIZE = 16*1024*1024
   integer, parameter :: CHUNK_OVERLAP = 256

   !> Local per-frame read window: starts small, grows (up to the cap) only
   !> for unusually large structures whose frame doesn't fit the first try
   integer, parameter :: FRAME_WINDOW_INIT = 65536
   integer, parameter :: FRAME_WINDOW_MAX  = 64*1024*1024

   character(len=*), parameter :: STRUCT_MARKER = '"structure": {"@module"'

   !> Per-config metadata stored after the pre-scan
   type, public :: json_config_info
      integer            :: natoms       = 0
      character(len=2)   :: class_id     = '01'   ! single class: MPtrj has no class=/subset= equivalent
      character(len=6)   :: cnumber      = '      '
      logical            :: has_energy   = .false.
      logical            :: has_forces   = .false.
      logical            :: has_stress   = .false.
      real(kind_double)  :: energy       = 0.d0   ! corrected_total_energy
      real(kind_double)  :: stress(9)    = 0.d0   ! row-major 3x3, as given in the file
      real(kind_double)  :: lattice(3,3) = 0.d0   ! columns = a,b,c
      integer(8)         :: struct_pos   = 0      ! stream pos (1-based) of structure's opening '{'
      integer(8)         :: frame_end    = 0      ! stream pos (1-based) of the frame's closing '}'
      logical            :: keep         = .true. ! .false. if dropped by drop_short_dist filtering
   end type json_config_info

   type(json_config_info), dimension(:), allocatable, public :: json_configs
   integer, public :: json_n_configs_total = 0
   integer, dimension(:), allocatable, public :: json_config_map

   !> Fixed-size, statically-allocated read buffers used for every stream
   !> positioned read in this module. This is deliberate: reading into an
   !> allocatable/deferred-length CHARACTER variable via a positioned stream
   !> read (the whole variable OR any slice of it) has been found to corrupt
   !> data - and sometimes the heap - with this gfortran version. Fixed-size
   !> buffers read correctly in all cases tested, so all I/O targets a slice
   !> of one of these two, and only plain (non-I/O) string assignment is
   !> ever used to hand data back through an allocatable CHARACTER result.
   !> Safe because all JSON reading in this module is rank0/subrank0-only
   !> and strictly sequential (no concurrent use of these buffers).
   character(len=CHUNK_SIZE), save :: chunk_buf
   character(len=FRAME_WINDOW_MAX), save :: frame_buf

   public :: json_detect_mode
   public :: json_scan_database
   public :: json_prepare_database
   public :: read_json_config

contains

   !---------------------------------------------------------------------------
   !> Detect whether db_path points to a JSON file
   !---------------------------------------------------------------------------
   subroutine json_detect_mode()
      use module_db_setup, only: db_path
      implicit none
      integer :: l

      db_json = .false.
      if (.not. allocated(db_path)) return
      l = len_trim(db_path)
      if (l < 5) return
      if (db_path(l-4:l) == '.json') db_json = .true.
   end subroutine json_detect_mode


   !---------------------------------------------------------------------------
   !> From buf(open_pos:open_pos) == open_char, find the matching close_char,
   !> honouring nesting depth and skipping the contents of "..." strings.
   !> Returns 0 if not found within buf.
   !---------------------------------------------------------------------------
   function find_matching(buf, open_pos, open_char, close_char) result(close_pos)
      implicit none
      character(len=*), intent(in) :: buf
      integer, intent(in) :: open_pos
      character(len=1), intent(in) :: open_char, close_char
      integer :: close_pos
      integer :: depth, i, n
      logical :: in_str, esc
      character(len=1) :: c

      close_pos = 0
      n = len(buf)
      depth = 1
      in_str = .false.
      esc = .false.
      i = open_pos + 1
      do while (i <= n)
         c = buf(i:i)
         if (in_str) then
            if (esc) then
               esc = .false.
            else if (c == '\') then
               esc = .true.
            else if (c == '"') then
               in_str = .false.
            end if
         else
            if (c == '"') then
               in_str = .true.
            else if (c == open_char) then
               depth = depth + 1
            else if (c == close_char) then
               depth = depth - 1
               if (depth == 0) then
                  close_pos = i
                  return
               end if
            end if
         end if
         i = i + 1
      end do
   end function find_matching


   !---------------------------------------------------------------------------
   !> Scan forward from from_pos (string-aware) for the first '}' - used to
   !> find a frame's own closing brace once its "structure" sub-object has
   !> already been consumed: the remaining flat fields (energy scalars,
   !> force/stress arrays, magmom, bandgap, mp_id string) contain no nested
   !> '{', only [...] arrays and strings, so the first unescaped '}' found
   !> here is unambiguously the frame's close. Returns 0 if not found.
   !---------------------------------------------------------------------------
   function find_frame_close(buf, from_pos) result(close_pos)
      implicit none
      character(len=*), intent(in) :: buf
      integer, intent(in) :: from_pos
      integer :: close_pos
      integer :: i, n
      logical :: in_str, esc
      character(len=1) :: c

      close_pos = 0
      n = len(buf)
      in_str = .false.
      esc = .false.
      i = from_pos
      do while (i <= n)
         c = buf(i:i)
         if (in_str) then
            if (esc) then
               esc = .false.
            else if (c == '\') then
               esc = .true.
            else if (c == '"') then
               in_str = .false.
            end if
         else
            if (c == '"') then
               in_str = .true.
            else if (c == '}') then
               close_pos = i
               return
            end if
         end if
         i = i + 1
      end do
   end function find_frame_close


   !---------------------------------------------------------------------------
   !> Locate a top-level "key": inside buf(from_pos:to_pos) (exact key match,
   !> since the full quoted key + colon is matched literally - "energy":
   !> will not match inside "corrected_total_energy": or "energy_per_atom":).
   !> Returns 0 if not found.
   !---------------------------------------------------------------------------
   function find_key(buf, key, from_pos, to_pos) result(kpos)
      implicit none
      character(len=*), intent(in) :: buf, key
      integer, intent(in) :: from_pos, to_pos
      integer :: kpos
      character(len=:), allocatable :: pattern
      integer :: i, plen

      pattern = '"'//trim(key)//'":'
      plen = len(pattern)
      kpos = 0
      if (to_pos - plen + 1 < from_pos) return
      do i = from_pos, to_pos - plen + 1
         if (buf(i:i+plen-1) == pattern) then
            kpos = i
            return
         end if
      end do
   end function find_key


   !---------------------------------------------------------------------------
   !> Extract a flat list of n reals from buf(start_pos:end_pos), which is
   !> expected to be a JSON number, or a (possibly nested) array of numbers
   !> such as [1,2,3] or [[1,2,3],[4,5,6]] - '[', ']' and ',' are blanked out
   !> and the remaining whitespace-separated numbers are read in file order.
   !---------------------------------------------------------------------------
   subroutine extract_flat_reals(buf, start_pos, end_pos, arr, n, ok)
      implicit none
      character(len=*), intent(in) :: buf
      integer, intent(in) :: start_pos, end_pos, n
      real(kind_double), intent(out) :: arr(n)
      logical, intent(out) :: ok
      character(len=:), allocatable :: tmp
      integer :: i, ios

      ok = .false.
      if (end_pos < start_pos) return
      tmp = buf(start_pos:end_pos)
      do i = 1, len(tmp)
         if (tmp(i:i) == '[' .or. tmp(i:i) == ']' .or. tmp(i:i) == ',') tmp(i:i) = ' '
      end do
      read(tmp, *, iostat=ios) arr(1:n)
      if (ios == 0) ok = .true.
   end subroutine extract_flat_reals


   !---------------------------------------------------------------------------
   !> Read a window of the file starting at 1-based stream position start_pos,
   !> growing it (up to FRAME_WINDOW_MAX) until buf holds a complete frame:
   !> the structure sub-object plus the flat fields through the frame's own
   !> closing '}'. Returns local positions (within buf) for the structure's
   !> open/close and the frame's close, or ok=.false. if parsing failed.
   !---------------------------------------------------------------------------
   subroutine read_one_frame(path, start_pos, buf, struct_open, struct_close, frame_close, ok)
      implicit none
      character(len=*), intent(in) :: path
      integer(8), intent(in) :: start_pos
      character(len=:), allocatable, intent(out) :: buf
      integer, intent(out) :: struct_open, struct_close, frame_close
      logical, intent(out) :: ok

      integer :: window_size, inp, ios

      ok = .false.
      struct_open = 0; struct_close = 0; frame_close = 0
      window_size = FRAME_WINDOW_INIT

      do
         frame_buf(1:window_size) = ' '
         open (newunit=inp, file=trim(path), status='old', action='read', &
               access='stream', form='unformatted', iostat=ios)
         if (ios /= 0) return
         read (inp, pos=start_pos, iostat=ios) frame_buf(1:window_size)
         close (inp)
         ! iostat may be nonzero at true EOF but frame_buf is still usable up
         ! to however many bytes remained; only bail out if it stayed blank
         struct_open = index(frame_buf(1:window_size), '{')
         if (struct_open > 0) then
            struct_close = find_matching(frame_buf(1:window_size), struct_open, '{', '}')
            if (struct_close > 0) then
               frame_close = find_frame_close(frame_buf(1:window_size), struct_close + 1)
               if (frame_close > 0) then
                  buf = frame_buf(1:window_size)
                  ok = .true.
                  return
               end if
            end if
         end if
         if (window_size >= FRAME_WINDOW_MAX) return
         window_size = min(window_size*4, FRAME_WINDOW_MAX)
      end do
   end subroutine read_one_frame


   !---------------------------------------------------------------------------
   !> Extract lattice (3x3, columns=a,b,c) from the structure sub-object
   !> buf(struct_open:struct_close). ok=.false. if no matrix key found.
   !---------------------------------------------------------------------------
   subroutine extract_lattice(buf, struct_open, struct_close, lattice, ok)
      implicit none
      character(len=*), intent(in) :: buf
      integer, intent(in) :: struct_open, struct_close
      real(kind_double), intent(out) :: lattice(3, 3)
      logical, intent(out) :: ok
      integer :: kpos, arr_open, arr_close
      real(kind_double) :: flat(9)

      lattice = 0.d0
      ok = .false.
      kpos = find_key(buf, 'matrix', struct_open, struct_close)
      if (kpos == 0) return
      arr_open = index(buf(kpos:struct_close), '[')
      if (arr_open == 0) return
      arr_open = kpos + arr_open - 1
      arr_close = find_matching(buf, arr_open, '[', ']')
      if (arr_close == 0) return
      call extract_flat_reals(buf, arr_open, arr_close, flat, 9, ok)
      if (.not. ok) return
      ! pymatgen "matrix" rows are the lattice vectors a,b,c (row-major flat
      ! order a1,a2,a3,b1,b2,b3,c1,c2,c3) -> store as columns to match the
      ! rest of Milady's convention (cell(:,1)=a, cell(:,2)=b, cell(:,3)=c)
      lattice(:, 1) = flat(1:3)
      lattice(:, 2) = flat(4:6)
      lattice(:, 3) = flat(7:9)
   end subroutine extract_lattice


   !---------------------------------------------------------------------------
   !> Extract, in site order, the element symbol and Cartesian xyz for every
   !> site inside the structure sub-object buf(struct_open:struct_close).
   !> Assumes single-occupancy sites (species list of length 1), which is
   !> the overwhelming common case in MPtrj; a site with anything else is
   !> reported via ok=.false. so the caller can skip/flag that config.
   !---------------------------------------------------------------------------
   subroutine extract_sites(buf, struct_open, struct_close, natoms, species, pos_cart, ok)
      implicit none
      character(len=*), intent(in) :: buf
      integer, intent(in) :: struct_open, struct_close, natoms
      character(len=2), intent(out) :: species(natoms)
      real(kind_double), intent(out) :: pos_cart(3, natoms)
      logical, intent(out) :: ok
      integer :: search_from, elem_kpos, elem_qstart, elem_qend
      integer :: xyz_kpos, xyz_open, xyz_close
      integer :: ia
      real(kind_double) :: xyz3(3)
      logical :: got

      ok = .false.
      search_from = struct_open
      do ia = 1, natoms
         elem_kpos = find_key(buf, 'element', search_from, struct_close)
         if (elem_kpos == 0) return
         elem_qstart = index(buf(elem_kpos:struct_close), '"', back=.false.)
         ! find_key already located the exact "element": token; the value's
         ! opening quote is the next '"' after the colon
         elem_qstart = elem_kpos + len('"element":')
         elem_qstart = index(buf(elem_qstart:struct_close), '"') + elem_qstart - 1
         if (elem_qstart == 0) return
         elem_qend = index(buf(elem_qstart+1:struct_close), '"')
         if (elem_qend == 0) return
         elem_qend = elem_qstart + elem_qend
         species(ia) = buf(elem_qstart+1:elem_qend-1)

         xyz_kpos = find_key(buf, 'xyz', elem_qend, struct_close)
         if (xyz_kpos == 0) return
         xyz_open = index(buf(xyz_kpos:struct_close), '[')
         if (xyz_open == 0) return
         xyz_open = xyz_kpos + xyz_open - 1
         xyz_close = find_matching(buf, xyz_open, '[', ']')
         if (xyz_close == 0) return
         call extract_flat_reals(buf, xyz_open, xyz_close, xyz3, 3, got)
         if (.not. got) return
         pos_cart(:, ia) = xyz3(:)

         search_from = xyz_close + 1
      end do
      ok = .true.
   end subroutine extract_sites


   !---------------------------------------------------------------------------
   !> Count sites (via "element": occurrences) inside the structure
   !> sub-object buf(struct_open:struct_close).
   !---------------------------------------------------------------------------
   function count_sites(buf, struct_open, struct_close) result(n)
      implicit none
      character(len=*), intent(in) :: buf
      integer, intent(in) :: struct_open, struct_close
      integer :: n
      integer :: kpos, from_pos

      n = 0
      from_pos = struct_open
      do
         kpos = find_key(buf, 'element', from_pos, struct_close)
         if (kpos == 0) exit
         n = n + 1
         from_pos = kpos + 1
      end do
   end function count_sites


   !---------------------------------------------------------------------------
   !> Pre-scan the JSON file: find every "structure" marker (rank 0 only,
   !> chunked overlapping search over the whole - possibly many-GB - file),
   !> then extract each frame's lightweight metadata (natoms, energy,
   !> lattice, stress, has_* flags, byte offsets) via a small local read.
   !> When drop_short_dist>0, atom positions are additionally parsed here to
   !> run the same short-pair filter used by the POSCAR/XYZ paths. Results
   !> are broadcast to all ranks via MPI_BCAST on mpi_comm_mld.
   !---------------------------------------------------------------------------
   subroutine json_scan_database()
      use module_db_setup, only: db_path, drop_short_dist
      use mld_mpi, only: mpi_comm_mld, rangml, mld_ierror
      use mpi
      use mld_logger
      implicit none

      integer(8) :: n_markers
      integer(8), allocatable :: marker_pos(:)
      integer :: n_total, i, n_dropped, n_bad
      logical :: ok

      _NAMECURRENT_("json_scan_database")
      _MLD_BEGIN_

      n_total = 0
      if (rangml == 0) then
         call json_find_markers(db_path, marker_pos, n_markers)
         write (6, '("ML: JSON database scan: found ", i0, " configurations in ", a)') &
            n_markers, trim(db_path)
      end if

      call MPI_BCAST(n_markers, 1, MPI_INTEGER8, 0, mpi_comm_mld, mld_ierror)
      n_total = int(n_markers)
      json_n_configs_total = n_total

      if (n_total == 0) then
         if (rangml == 0) write (6, *) 'ML: FATAL: no configurations found in JSON file'
         call MPI_Abort(mpi_comm_mld, 1, mld_ierror)
      end if

      if (allocated(json_configs)) deallocate (json_configs)
      allocate (json_configs(n_total))

      if (rangml == 0) then
         n_dropped = 0
         n_bad = 0
         do i = 1, n_total
            call json_scan_one(db_path, marker_pos(i), drop_short_dist, json_configs(i), ok)
            if (.not. ok) then
               n_bad = n_bad + 1
               json_configs(i)%keep = .false.
               cycle
            end if
            if (.not. json_configs(i)%keep) n_dropped = n_dropped + 1
            write (json_configs(i)%cnumber, '(i6.6)') i
         end do
         if (n_bad > 0) then
            write (6, '("ML: WARNING: ", i0, " of ", i0, &
               &" JSON configurations could not be parsed and were excluded")') n_bad, n_total
         end if
         if (drop_short_dist > 0.d0) then
            write (6, '("ML: drop_short_dist = ", f10.5, " Angstrom: dropped ", i0, &
               &" of ", i0, " JSON configs (interatomic distance below threshold)")') &
               drop_short_dist, n_dropped, n_total
         end if
         deallocate (marker_pos)
      end if

      call json_bcast_all_configs(n_total)

      _MLD_END_
   end subroutine json_scan_database


   !---------------------------------------------------------------------------
   !> Bulk streaming pass: find the byte offset of every "structure": marker
   !> in the file, using fixed-size overlapping chunk reads so a marker split
   !> across a chunk boundary is never missed and never double-counted.
   !---------------------------------------------------------------------------
   subroutine json_find_markers(path, marker_pos, n_markers)
      implicit none
      character(len=*), intent(in) :: path
      integer(8), allocatable, intent(out) :: marker_pos(:)
      integer(8), intent(out) :: n_markers

      integer :: inp, ios
      integer(8) :: file_size, offset, read_len, cap
      integer :: mlen
      logical :: last_chunk

      mlen = len(STRUCT_MARKER)
      cap = 1024
      allocate (marker_pos(cap))
      n_markers = 0

      inquire (file=trim(path), size=file_size)
      if (file_size <= 0) return

      open (newunit=inp, file=trim(path), status='old', action='read', &
            access='stream', form='unformatted', iostat=ios)
      if (ios /= 0) return

      offset = 0
      do
         read_len = min(int(CHUNK_SIZE, 8), file_size - offset)
         if (read_len <= 0) exit
         last_chunk = (offset + read_len >= file_size)
         chunk_buf(1:read_len) = ' '
         read (inp, pos=offset+1, iostat=ios) chunk_buf(1:read_len)
         if (ios /= 0 .and. .not. last_chunk) exit

         call scan_chunk_markers(chunk_buf, int(read_len), mlen, &
            merge(int(read_len), int(read_len - CHUNK_OVERLAP), last_chunk), &
            offset, marker_pos, n_markers, cap)

         if (last_chunk) exit
         offset = offset + read_len - CHUNK_OVERLAP
      end do
      close (inp)
   end subroutine json_find_markers


   !---------------------------------------------------------------------------
   !> Find every occurrence of STRUCT_MARKER within buf(1:read_len), but only
   !> accept matches starting at local position <= limit (the overlap tail is
   !> left for the next chunk to pick up, avoiding double counting). Appends
   !> absolute (0-based) byte offsets to marker_pos, growing it as needed.
   !---------------------------------------------------------------------------
   subroutine scan_chunk_markers(buf, read_len, mlen, limit, offset, marker_pos, n_markers, cap)
      implicit none
      character(len=*), intent(in) :: buf
      integer, intent(in) :: read_len, mlen, limit
      integer(8), intent(in) :: offset
      integer(8), allocatable, intent(inout) :: marker_pos(:)
      integer(8), intent(inout) :: n_markers
      integer(8), intent(inout) :: cap
      integer :: search_from, found
      integer(8), allocatable :: tmp(:)

      search_from = 1
      do
         if (search_from > read_len) exit
         found = index(buf(search_from:read_len), STRUCT_MARKER)
         if (found == 0) exit
         found = found + search_from - 1
         if (found > limit) exit
         n_markers = n_markers + 1
         if (n_markers > cap) then
            allocate (tmp(cap))
            tmp = marker_pos(1:cap)
            cap = cap*2
            deallocate (marker_pos)
            allocate (marker_pos(cap))
            marker_pos(1:n_markers-1) = tmp(1:n_markers-1)
            deallocate (tmp)
         end if
         marker_pos(n_markers) = offset + int(found, 8) - 1_8   ! 0-based
         search_from = found + mlen
      end do
   end subroutine scan_chunk_markers


   !---------------------------------------------------------------------------
   !> Extract one configuration's scan-time metadata given the 0-based byte
   !> offset of its "structure": marker.
   !---------------------------------------------------------------------------
   subroutine json_scan_one(path, marker_abs0, drop_short_dist, cfg, ok)
      implicit none
      character(len=*), intent(in) :: path
      integer(8), intent(in) :: marker_abs0
      real(kind_double), intent(in) :: drop_short_dist
      type(json_config_info), intent(inout) :: cfg
      logical, intent(out) :: ok

      character(len=:), allocatable :: buf
      integer :: struct_open, struct_close, frame_close
      integer(8) :: start_pos
      logical :: got
      integer :: kpos, val_end
      character(len=2), allocatable :: species(:)
      real(kind_double), allocatable :: pos_cart(:, :)
      logical :: too_close

      ok = .false.
      start_pos = marker_abs0 + 1  ! 1-based stream pos
      call read_one_frame(path, start_pos, buf, struct_open, struct_close, frame_close, got)
      if (.not. got) return

      cfg%struct_pos = start_pos + int(struct_open, 8) - 1_8
      cfg%frame_end = start_pos + int(frame_close, 8) - 1_8

      call extract_lattice(buf, struct_open, struct_close, cfg%lattice, got)
      if (.not. got) return

      cfg%natoms = count_sites(buf, struct_open, struct_close)
      if (cfg%natoms <= 0) return

      ! energy target: corrected_total_energy
      kpos = find_key(buf, 'corrected_total_energy', struct_close, frame_close)
      if (kpos > 0) then
         call extract_scalar_real(buf, kpos + len('"corrected_total_energy":'), frame_close, &
                                    cfg%energy, cfg%has_energy)
      end if

      ! force presence (value is either an array or the literal null)
      kpos = find_key(buf, 'force', struct_close, frame_close)
      if (kpos > 0) then
         val_end = kpos + len('"force":')
         cfg%has_forces = (index(buf(val_end:min(val_end+10, frame_close)), 'null') == 0)
      end if

      ! stress presence + value
      kpos = find_key(buf, 'stress', struct_close, frame_close)
      if (kpos > 0) then
         val_end = kpos + len('"stress":')
         if (index(buf(val_end:min(val_end+10, frame_close)), 'null') == 0) then
            block
               integer :: arr_open, arr_close
               arr_open = index(buf(val_end:frame_close), '[')
               if (arr_open > 0) then
                  arr_open = val_end + arr_open - 1
                  arr_close = find_matching(buf, arr_open, '[', ']')
                  if (arr_close > 0) then
                     call extract_flat_reals(buf, arr_open, arr_close, cfg%stress, 9, got)
                     cfg%has_stress = got
                  end if
               end if
            end block
         end if
      end if

      if (drop_short_dist > 0.d0) then
         allocate (species(cfg%natoms), pos_cart(3, cfg%natoms))
         call extract_sites(buf, struct_open, struct_close, cfg%natoms, species, pos_cart, got)
         if (got) then
            call config_has_short_pair(cfg%natoms, pos_cart, cfg%lattice, drop_short_dist, too_close)
            cfg%keep = .not. too_close
         end if
         deallocate (species, pos_cart)
      end if

      ok = .true.
   end subroutine json_scan_one


   !---------------------------------------------------------------------------
   !> Extract a bare numeric value starting at from_pos (first non-space
   !> character), stopping at the first ',' or '}' at bracket depth 0.
   !---------------------------------------------------------------------------
   subroutine extract_scalar_real(buf, from_pos, to_pos, val, found)
      implicit none
      character(len=*), intent(in) :: buf
      integer, intent(in) :: from_pos, to_pos
      real(kind_double), intent(out) :: val
      logical, intent(out) :: found
      integer :: i, vend, ios

      found = .false.
      val = 0.d0
      i = from_pos
      do while (i <= to_pos .and. buf(i:i) == ' ')
         i = i + 1
      end do
      if (i > to_pos) return
      if (buf(i:min(i+3, to_pos)) == 'null') return
      vend = i
      do while (vend <= to_pos)
         if (buf(vend:vend) == ',' .or. buf(vend:vend) == '}') exit
         vend = vend + 1
      end do
      read (buf(i:vend-1), *, iostat=ios) val
      if (ios == 0) found = .true.
   end subroutine extract_scalar_real


   !---------------------------------------------------------------------------
   !> Broadcast json_configs(1:n) from rank 0 to all ranks.
   !---------------------------------------------------------------------------
   subroutine json_bcast_all_configs(n)
      use mld_mpi, only: mpi_comm_mld, mld_ierror
      use mpi
      implicit none
      integer, intent(in) :: n
      integer :: i, idx
      integer, allocatable :: ibuf(:), lbuf(:)
      integer(8), allocatable :: i8buf(:)
      real(kind_double), allocatable :: rbuf(:)

      if (n == 0) return

      allocate (ibuf(n))
      ibuf = json_configs(:)%natoms
      call MPI_BCAST(ibuf, n, MPI_INTEGER, 0, mpi_comm_mld, mld_ierror)
      json_configs(:)%natoms = ibuf
      deallocate (ibuf)

      allocate (i8buf(2*n))
      do i = 1, n
         i8buf(i) = json_configs(i)%struct_pos
         i8buf(n+i) = json_configs(i)%frame_end
      end do
      call MPI_BCAST(i8buf, 2*n, MPI_INTEGER8, 0, mpi_comm_mld, mld_ierror)
      do i = 1, n
         json_configs(i)%struct_pos = i8buf(i)
         json_configs(i)%frame_end = i8buf(n+i)
      end do
      deallocate (i8buf)

      allocate (lbuf(4*n))
      do i = 1, n
         lbuf(i)     = merge(1, 0, json_configs(i)%has_energy)
         lbuf(n+i)   = merge(1, 0, json_configs(i)%has_forces)
         lbuf(2*n+i) = merge(1, 0, json_configs(i)%has_stress)
         lbuf(3*n+i) = merge(1, 0, json_configs(i)%keep)
      end do
      call MPI_BCAST(lbuf, 4*n, MPI_INTEGER, 0, mpi_comm_mld, mld_ierror)
      do i = 1, n
         json_configs(i)%has_energy  = (lbuf(i)     == 1)
         json_configs(i)%has_forces  = (lbuf(n+i)   == 1)
         json_configs(i)%has_stress  = (lbuf(2*n+i) == 1)
         json_configs(i)%keep        = (lbuf(3*n+i) == 1)
      end do
      deallocate (lbuf)

      allocate (rbuf(19*n))
      do i = 1, n
         idx = (i-1)*19
         rbuf(idx+1) = json_configs(i)%energy
         rbuf(idx+2:idx+10) = json_configs(i)%stress(1:9)
         rbuf(idx+11) = json_configs(i)%lattice(1,1)
         rbuf(idx+12) = json_configs(i)%lattice(2,1)
         rbuf(idx+13) = json_configs(i)%lattice(3,1)
         rbuf(idx+14) = json_configs(i)%lattice(1,2)
         rbuf(idx+15) = json_configs(i)%lattice(2,2)
         rbuf(idx+16) = json_configs(i)%lattice(3,2)
         rbuf(idx+17) = json_configs(i)%lattice(1,3)
         rbuf(idx+18) = json_configs(i)%lattice(2,3)
         rbuf(idx+19) = json_configs(i)%lattice(3,3)
      end do
      call MPI_BCAST(rbuf, 19*n, MPI_DOUBLE_PRECISION, 0, mpi_comm_mld, mld_ierror)
      do i = 1, n
         idx = (i-1)*19
         json_configs(i)%energy = rbuf(idx+1)
         json_configs(i)%stress(1:9) = rbuf(idx+2:idx+10)
         json_configs(i)%lattice(1,1) = rbuf(idx+11)
         json_configs(i)%lattice(2,1) = rbuf(idx+12)
         json_configs(i)%lattice(3,1) = rbuf(idx+13)
         json_configs(i)%lattice(1,2) = rbuf(idx+14)
         json_configs(i)%lattice(2,2) = rbuf(idx+15)
         json_configs(i)%lattice(3,2) = rbuf(idx+16)
         json_configs(i)%lattice(1,3) = rbuf(idx+17)
         json_configs(i)%lattice(2,3) = rbuf(idx+18)
         json_configs(i)%lattice(3,3) = rbuf(idx+19)
      end do
      deallocate (rbuf)
   end subroutine json_bcast_all_configs


   !---------------------------------------------------------------------------
   !> Set up db_model/config_real/etc. from the JSON scan. Mirrors
   !> xyz_prepare_database, but with a single fixed class ('01').
   !---------------------------------------------------------------------------
   subroutine json_prepare_database()
      use ml_in_ndm_module, only: rangml, debug
      use module_db_setup, only: db_file, db_train, db_test, &
         iconf_data, iconf_data_train, iconf_data_test
      use derived_types, only: db_model, config_real, config_desc
      use mld_unit
      use mld_logger
      use mld_mpi
      implicit none

      integer :: i, n_total, n_kept, iunit, ierr, n_db_class
      logical :: lexist

      _NAMECURRENT_("json_prepare_database")
      _MLD_BEGIN_

      n_total = json_n_configs_total
      n_kept = count(json_configs(1:n_total)%keep)

      if (rangml == 0) then
         write (6, '("ML: JSON class 01: ", i0, " configs in JSON file")') n_kept
      end if

      inquire (file=db_file, exist=lexist)
      if (.not. lexist) then
         call log_critical('json_prepare_database: db_file not found: '//trim(db_file))
         call mld_mpi_finalize('json_prepare_database: provide db_model.in')
         stop
      end if
      open (newunit=iunit, file=db_file, status='old', action='read', iostat=ierr)
      if (ierr /= 0) then
         call log_critical('json_prepare_database: error opening db_file: '//trim(db_file))
         call mld_mpi_finalize('json_prepare_database: check file')
         stop
      end if
      call read_n_db_class(iunit, n_db_class)
      call read_db_file(iunit, n_db_class)
      close (iunit)

      do i = 1, size(db_model)
         if (db_model(i)%class /= '01') then
            if (rangml == 0) write (6, '("ML: ERROR: class ", a, &
               &" in db_model.in not found in JSON file (single class 01 only)")') db_model(i)%class
            call mld_mpi_finalize('json_prepare_database: class not found in JSON')
            stop
         end if
         if (db_model(i)%no_total > n_kept) then
            if (rangml == 0) then
               write (6, '("ML: WARNING: class 01 db_model no_total=", i0, &
                  &" exceeds JSON count=", i0, ", clamping")') db_model(i)%no_total, n_kept
            end if
            db_model(i)%no_total = n_kept
            if (db_model(i)%no_selec > db_model(i)%no_total) then
               db_model(i)%no_selec = db_model(i)%no_total
            end if
         end if
      end do

      iconf_data_train = 0
      iconf_data = 0
      do i = 1, size(db_model)
         iconf_data_train = iconf_data_train + db_model(i)%no_selec
         iconf_data = iconf_data + db_model(i)%no_total - db_model(i)%no_start + 1
      end do
      iconf_data_test = iconf_data - iconf_data_train

      if (allocated(config_real)) deallocate (config_real); allocate (config_real(iconf_data))
      if (allocated(config_desc)) deallocate (config_desc); allocate (config_desc(iconf_data))
      if (allocated(db_train)) deallocate (db_train); allocate (db_train(iconf_data_train))
      if (allocated(db_test)) deallocate (db_test); allocate (db_test(max(iconf_data_test, 0)))

      call prepare_name_file_from_db()

      ! config_real(i)%no_file_in_db_line is the 1-based sequence number
      ! among KEPT json_configs entries; build the map directly from that
      if (allocated(json_config_map)) deallocate (json_config_map)
      allocate (json_config_map(iconf_data))
      block
         integer :: j, kept_seen
         kept_seen = 0
         do i = 1, iconf_data
            do while (kept_seen < n_total)
               kept_seen = kept_seen + 1
               if (json_configs(kept_seen)%keep) exit
            end do
            json_config_map(i) = kept_seen
         end do
      end block

      do i = 1, iconf_data
         config_real(i)%has_energy = config_real(i)%has_energy .and. json_configs(json_config_map(i))%has_energy
         config_real(i)%has_force  = config_real(i)%has_force  .and. json_configs(json_config_map(i))%has_forces
         config_real(i)%has_stress = config_real(i)%has_stress .and. json_configs(json_config_map(i))%has_stress
      end do

      call get_distinct_classes()

      if (rangml == 0) then
         write (6, '("ML: JSON database prepared: ", i0, " total, ", i0, " train, ", i0, " test")') &
            iconf_data, iconf_data_train, iconf_data_test
      end if

      _MLD_END_
   end subroutine json_prepare_database


   !---------------------------------------------------------------------------
   !> Read one configuration's atoms (species, Cartesian positions, forces)
   !> from the JSON file, using the byte range recorded during the scan.
   !> Only subrank==0 reads; results are broadcast to the subworld.
   !---------------------------------------------------------------------------
   subroutine read_json_config(ifile)
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
      integer :: ijson, natoms, i, ij, ii, i_p, icnt, im_local, nb_elements
      character(len=:), allocatable :: buf
      integer :: struct_open, struct_close, frame_close
      logical :: got, ok_read
      real(kind_double) :: box(3,3), box_inv(3,3), st(6)
      character(len=2), allocatable :: atom_species(:)
      real(kind_double), allocatable :: atom_pos(:,:), atom_force(:,:)
      character(len=2), allocatable :: element_list(:)
      integer, allocatable :: nspecies(:), ityp(:)
      real(kind_double), allocatable :: xp(:,:), xc(:,:), fp(:,:), l_spin(:,:)
      real(kind_double), allocatable :: vect_buffer(:)
      character(len=2) :: ch_buffer
      integer :: n_unique, found
      integer :: kpos, val_end, arr_open, arr_close
      real(kind_double) :: E_total, volume

      _NAMECURRENT_("read_json_config")
      _MLD_BEGIN_

      ijson = json_config_map(ifile)
      natoms = json_configs(ijson)%natoms
      box = json_configs(ijson)%lattice
      st(1) = json_configs(ijson)%stress(1)
      st(2) = json_configs(ijson)%stress(5)
      st(3) = json_configs(ijson)%stress(9)
      st(4) = json_configs(ijson)%stress(2)
      st(5) = json_configs(ijson)%stress(3)
      st(6) = json_configs(ijson)%stress(6)
      E_total = json_configs(ijson)%energy

      allocate (atom_species(natoms), atom_pos(3,natoms), atom_force(3,natoms))
      atom_species = '  '
      atom_pos = 0.d0
      atom_force = 0.d0

      if (subrank == 0) then
         block
            integer :: window
            integer :: inp, ios
            ! frame_end-struct_pos+1 <= FRAME_WINDOW_MAX always holds: the
            ! scan phase (read_one_frame) only recorded these offsets for a
            ! frame it successfully bounded within that same cap.
            window = int(json_configs(ijson)%frame_end - json_configs(ijson)%struct_pos + 1)
            frame_buf(1:window) = ' '
            open (newunit=inp, file=trim(db_path), status='old', action='read', &
                  access='stream', form='unformatted', iostat=ios)
            if (ios /= 0) then
               write (6,*) 'ML: FATAL: cannot open JSON file: ', trim(db_path)
               stop 'FATAL: cannot open JSON file in read_json_config'
            end if
            read (inp, pos=json_configs(ijson)%struct_pos, iostat=ios) frame_buf(1:window)
            close (inp)
            buf = frame_buf(1:window)
         end block
         struct_open = 1
         struct_close = find_matching(buf, struct_open, '{', '}')
         frame_close = len(buf)
         call extract_sites(buf, struct_open, struct_close, natoms, atom_species, atom_pos, ok_read)
         if (.not. ok_read) then
            write (6,*) 'ML: FATAL: could not re-parse sites for JSON config ', ifile
            stop 'FATAL: JSON site parse error in read_json_config'
         end if

         if (json_configs(ijson)%has_forces) then
            kpos = find_key(buf, 'force', struct_close, frame_close)
            if (kpos > 0) then
               val_end = kpos + len('"force":')
               arr_open = index(buf(val_end:frame_close), '[')
               if (arr_open > 0) then
                  arr_open = val_end + arr_open - 1
                  arr_close = find_matching(buf, arr_open, '[', ']')
                  if (arr_close > 0) then
                     block
                        real(kind_double) :: fflat(3*natoms)
                        call extract_flat_reals(buf, arr_open, arr_close, fflat, 3*natoms, got)
                        if (got) then
                           do i = 1, natoms
                              atom_force(:, i) = fflat(3*(i-1)+1:3*(i-1)+3)
                           end do
                        end if
                     end block
                  end if
               end if
            end if
         end if
      end if

      ! ---- Broadcast atom data to subworld ----
      do i = 1, natoms
         ch_buffer = atom_species(i)
         call my_broadcast_char(ch_buffer, 0, subworld, codeml)
         atom_species(i) = ch_buffer
      end do

      allocate (vect_buffer(natoms))
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
      deallocate (vect_buffer)

      allocate (vect_buffer(9))
      vect_buffer(1:3) = box(:,1); vect_buffer(4:6) = box(:,2); vect_buffer(7:9) = box(:,3)
      call my_broadcast_vect_real(vect_buffer, 0, subworld, codeml)
      box(:,1) = vect_buffer(1:3); box(:,2) = vect_buffer(4:6); box(:,3) = vect_buffer(7:9)
      deallocate (vect_buffer)

      allocate (vect_buffer(6))
      vect_buffer(1:6) = st(1:6)
      call my_broadcast_vect_real(vect_buffer, 0, subworld, codeml)
      st(1:6) = vect_buffer(1:6)
      deallocate (vect_buffer)

      call my_broadcast_val_real(E_total, 0, subworld, codeml)

      ! ---- Unique species, sorted, POSCAR-style grouping ----
      n_unique = 0
      allocate (element_list(natoms))
      element_list = '  '
      do i = 1, natoms
         found = 0
         do ij = 1, n_unique
            if (element_list(ij) == atom_species(i)) then
               found = ij; exit
            end if
         end do
         if (found == 0) then
            n_unique = n_unique + 1
            element_list(n_unique) = atom_species(i)
         end if
      end do
      call sort_species_json(element_list, n_unique)

      nb_elements = n_unique
      im_local = natoms

      allocate (nspecies(nb_elements))
      nspecies = 0
      do i = 1, natoms
         do ij = 1, nb_elements
            if (atom_species(i) == element_list(ij)) then
               nspecies(ij) = nspecies(ij) + 1
               exit
            end if
         end do
      end do

      allocate (ityp(im_local))
      allocate (xp(3,im_local), xc(3,im_local), fp(3,im_local), l_spin(3,im_local))
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

      if (abs(box(1,1)) > 1.d-10 .or. abs(box(2,2)) > 1.d-10 .or. abs(box(3,3)) > 1.d-10 .or. &
          abs(box(1,2)) > 1.d-10 .or. abs(box(2,1)) > 1.d-10 .or. abs(box(1,3)) > 1.d-10) then
         call matinv_gen(box, box_inv)
         xc(1:3,1:im_local) = matmul(box_inv(1:3,1:3), xp(1:3,1:im_local))
      else
         box = 0.d0
         box(1,1) = 100.d0; box(2,2) = 100.d0; box(3,3) = 100.d0
         call matinv_gen(box, box_inv)
         xc(1:3,1:im_local) = matmul(box_inv(1:3,1:3), xp(1:3,1:im_local))
      end if

      if (im_local > imm) then
         imm = im_local
         if (debug) write (6,'("ML: WARNING imm increased to ",i7," in read_json_config")') imm
      end if
      im = im_local

      if (allocated(config_real(ifile)%ref_energy_per_element)) deallocate (config_real(ifile)%ref_energy_per_element)
      allocate (config_real(ifile)%ref_energy_per_element(nb_elements))
      if (allocated(config_real(ifile)%mass_per_type)) deallocate (config_real(ifile)%mass_per_type)
      allocate (config_real(ifile)%mass_per_type(nb_elements))
      if (allocated(config_real(ifile)%Z_per_type)) deallocate (config_real(ifile)%Z_per_type)
      allocate (config_real(ifile)%Z_per_type(nb_elements))
      if (allocated(config_real(ifile)%covalent_radius_per_type)) deallocate (config_real(ifile)%covalent_radius_per_type)
      allocate (config_real(ifile)%covalent_radius_per_type(nb_elements))
      if (allocated(config_real(ifile)%fix_type_poscar_to_periodic)) deallocate (config_real(ifile)%fix_type_poscar_to_periodic)
      allocate (config_real(ifile)%fix_type_poscar_to_periodic(nb_elements))
      if (allocated(config_real(ifile)%itype_to_global)) deallocate (config_real(ifile)%itype_to_global)
      allocate (config_real(ifile)%itype_to_global(nb_elements))

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
               write (6,*) "ML: Unknown element in JSON database: ", element_list(ij)
               write (6,*) "ML: Please update chemical_elements in input file."
               stop 'fatal in read_json_config, unknown element (weighted)'
            end if
         end do
         if (linvisible) then
            if (allocated(config_real(ifile)%invisible_per_type)) deallocate (config_real(ifile)%invisible_per_type)
            allocate (config_real(ifile)%invisible_per_type(nb_elements))
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
               write (6,*) "ML: Unknown element in JSON database: ", element_list(ij)
               write (6,*) "ML: Please update chemical_elements."
               stop 'fatal in read_json_config, unknown element'
            end if
         end do
      end if

      if (weighted) then
         fix_weighted_for_element_ini(:) = fix_weighted_for_element(:)
         if (weighted_3ch) fix_weighted_for_element_3ch_ini(:) = fix_weighted_for_element_3ch(:)
      end if

      if (allocated(config_real(ifile)%itype)) deallocate (config_real(ifile)%itype)
      allocate (config_real(ifile)%itype(im_local))
      if (allocated(config_real(ifile)%itype_db)) deallocate (config_real(ifile)%itype_db)
      allocate (config_real(ifile)%itype_db(im_local))
      if (allocated(config_real(ifile)%pos_cart)) deallocate (config_real(ifile)%pos_cart)
      allocate (config_real(ifile)%pos_cart(3,im_local))
      if (allocated(config_real(ifile)%pos_crst)) deallocate (config_real(ifile)%pos_crst)
      allocate (config_real(ifile)%pos_crst(3,im_local))
      if (allocated(config_real(ifile)%force)) deallocate (config_real(ifile)%force)
      allocate (config_real(ifile)%force(3,im_local))
      if (allocated(config_real(ifile)%atomic_spin)) deallocate (config_real(ifile)%atomic_spin)
      allocate (config_real(ifile)%atomic_spin(3,im_local))

      config_real(ifile)%ntypes = nb_elements
      config_real(ifile)%nat = im_local
      config_real(ifile)%itype(1:im_local) = ityp(1:im_local)
      do ii = 1, im_local
         config_real(ifile)%itype_db(ii) = config_real(ifile)%itype_to_global(ityp(ii))
      end do
      config_real(ifile)%pos_crst(1:3,1:im_local) = xc(1:3,1:im_local)
      config_real(ifile)%pos_cart(1:3,1:im_local) = xp(1:3,1:im_local)
      config_real(ifile)%force(1:3,1:im_local) = fp(1:3,1:im_local)
      config_real(ifile)%atomic_spin = 0.d0
      config_real(ifile)%cell = box
      call calc_volume(config_real(ifile)%cell(:,1), config_real(ifile)%cell(:,2), &
                        config_real(ifile)%cell(:,3), volume)
      config_real(ifile)%volume = volume
      config_real(ifile)%energy(1) = E_total
      config_real(ifile)%energy(2) = E_total
      config_real(ifile)%energy(3) = E_total
      config_real(ifile)%stress(1:6) = st(1:6)

      _MLD_END_
   end subroutine read_json_config


   !> Sort species names alphabetically (bubble sort, n is small)
   subroutine sort_species_json(arr, n)
      implicit none
      character(len=2), intent(inout) :: arr(:)
      integer, intent(in) :: n
      integer :: i, j
      character(len=2) :: tmp
      logical :: swapped

      do i = 1, n - 1
         swapped = .false.
         do j = 1, n - i
            if (arr(j) > arr(j+1)) then
               tmp = arr(j); arr(j) = arr(j+1); arr(j+1) = tmp
               swapped = .true.
            end if
         end do
         if (.not. swapped) exit
      end do
   end subroutine sort_species_json

end module module_json
