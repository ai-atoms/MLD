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

!-------------------------------------------------------------------------------
! module_pair_cutoffs
!
! Setup and management of pair-specific inner cutoff and ZBL parameters.
!
! Logic:
!   - r_cut_in > 0  =>  global mode: all pairs get the same value (backward compatible)
!   - r_cut_in < 0  =>  automatic mode: computed from covalent radii per element pair
!   - Override strings (r_cut_in_pair_force, r1_zbl_pair_force, r2_zbl_pair_force)
!     can override specific pairs after the global/automatic initialization.
!
! Arrays populated:
!   r_cut_pair_in(i,j)  in module_neigh_local
!   r1_zbl_pair(i,j)    in module_zbl
!   r2_zbl_pair(i,j)    in module_zbl
!   rr_k2b_pair(i,j)    in module_zbl
!-------------------------------------------------------------------------------

module module_pair_cutoffs
   use module_kind_variables, only: kind_double
   implicit none

   ! override strings from namelist (format: "El1 El2 val ; El1 El2 val")
   character(len=5000) :: r_cut_in_pair_force  = ""
   character(len=5000) :: r1_zbl_pair_force    = ""
   character(len=5000) :: r2_zbl_pair_force    = ""

   ! automatic-mode parameters
   real(kind_double), parameter :: auto_scale_cov     = 0.60d0   ! fraction of sum of covalent radii
   real(kind_double), parameter :: auto_clamp_min     = 0.60d0   ! Ang, lower bound
   real(kind_double), parameter :: auto_clamp_max     = 1.80d0   ! Ang, upper bound
   real(kind_double), parameter :: auto_delta1_zbl    = 0.40d0   ! r1_zbl = r_cut_in - delta1
   real(kind_double), parameter :: auto_delta2_zbl    = 0.70d0   ! r2_zbl = r_cut_in + delta2
   real(kind_double), parameter :: auto_delta1_min    = 0.10d0   ! minimum r1_zbl

contains

   !----------------------------------------------------------------------------
   ! setup_pair_cutoffs
   !
   ! Must be called AFTER fix_type_of_atoms (species arrays are allocated).
   ! Populates r_cut_pair_in, r1_zbl_pair, r2_zbl_pair, rr_k2b_pair.
   !----------------------------------------------------------------------------
   subroutine setup_pair_cutoffs()
      use module_neigh_local, only: r_cut_in, r_cut_width_in, r_cut_pair_in
      use module_zbl, only: zbl_potential, zbl_type, zbl_mode_default_k2b, &
                            r1_zbl, r2_zbl, rr_k2b, &
                            r1_zbl_pair, r2_zbl_pair, rr_k2b_pair
      use module_chemical_species, only: fix_no_of_elements, fix_ch_elements, &
                                         fix_covalent_radius_elements
      use mld_logger
      use mld_string

      implicit none
      integer :: ii, jj, nsp
      real(kind_double) :: rcov_i, rcov_j, rval
      character(len=:), allocatable :: msg

      nsp = fix_no_of_elements

      ! --- allocate pair arrays ---
      if (allocated(r_cut_pair_in)) deallocate(r_cut_pair_in)
      allocate(r_cut_pair_in(nsp, nsp))

      if (allocated(r1_zbl_pair)) deallocate(r1_zbl_pair)
      allocate(r1_zbl_pair(nsp, nsp))

      if (allocated(r2_zbl_pair)) deallocate(r2_zbl_pair)
      allocate(r2_zbl_pair(nsp, nsp))

      if (allocated(rr_k2b_pair)) deallocate(rr_k2b_pair)
      allocate(rr_k2b_pair(nsp, nsp))

      ! =====================================================================
      ! 1. Initialize r_cut_pair_in
      ! =====================================================================
      if (r_cut_in > 0.d0) then
         ! GLOBAL MODE: same value for all pairs
         r_cut_pair_in(:,:) = r_cut_in
         call log_info('Pair cutoffs: GLOBAL mode for r_cut_in = '//vtoa(r_cut_in)//' Ang (same for all element pairs)')
      else
         ! AUTOMATIC MODE: element-specific from covalent radii
         call log_info('Pair cutoffs: AUTOMATIC mode for r_cut_in (computed from covalent radii)')
         do ii = 1, nsp
            rcov_i = fix_covalent_radius_elements(ii)  ! in pm
            do jj = 1, nsp
               rcov_j = fix_covalent_radius_elements(jj)  ! in pm
               rval = auto_scale_cov * (rcov_i + rcov_j) / 100.d0   ! pm -> Ang
               ! clamp
               rval = max(auto_clamp_min, min(auto_clamp_max, rval))
               r_cut_pair_in(ii, jj) = rval
            end do
         end do
      end if

      ! =====================================================================
      ! 2. Apply overrides from r_cut_in_pair_force
      ! =====================================================================
      if (len_trim(r_cut_in_pair_force) > 0) then
         call apply_pair_override(r_cut_in_pair_force, r_cut_pair_in, nsp, 'r_cut_in_pair_force')
      end if

      ! =====================================================================
      ! 3. Initialize r1_zbl_pair, r2_zbl_pair
      ! =====================================================================
      if (zbl_potential) then
         if (r1_zbl > 0.d0 .and. r_cut_in > 0.d0) then
            ! GLOBAL MODE: same value for all pairs
            r1_zbl_pair(:,:) = r1_zbl
            r2_zbl_pair(:,:) = r2_zbl
            call log_info('Pair cutoffs: GLOBAL mode for r1_zbl = '//vtoa(r1_zbl)//' r2_zbl = '//vtoa(r2_zbl))
         else
            ! AUTOMATIC MODE: derived from r_cut_pair_in
            call log_info('Pair cutoffs: AUTOMATIC mode for r1_zbl, r2_zbl (derived from pair-specific r_cut_in)')
            do ii = 1, nsp
               do jj = 1, nsp
                  rval = r_cut_pair_in(ii, jj)
                  r1_zbl_pair(ii, jj) = max(auto_delta1_min, rval - auto_delta1_zbl)
                  r2_zbl_pair(ii, jj) = rval + auto_delta2_zbl
               end do
            end do
         end if

         ! Apply overrides from r1_zbl_pair_force, r2_zbl_pair_force
         if (len_trim(r1_zbl_pair_force) > 0) then
            call apply_pair_override(r1_zbl_pair_force, r1_zbl_pair, nsp, 'r1_zbl_pair_force')
         end if
         if (len_trim(r2_zbl_pair_force) > 0) then
            call apply_pair_override(r2_zbl_pair_force, r2_zbl_pair, nsp, 'r2_zbl_pair_force')
         end if

         ! set rr_k2b_pair = r_cut_pair_in (transition point for bridge)
         rr_k2b_pair(:,:) = r_cut_pair_in(:,:)

         ! Update global scalars for backward compatibility:
         ! use minimum across all pairs so neighbor lists and validation still work
         r_cut_in = minval(r_cut_pair_in)
         r1_zbl   = minval(r1_zbl_pair)
         r2_zbl   = maxval(r2_zbl_pair)
         rr_k2b   = r_cut_in
      else
         ! No ZBL: r1/r2_zbl_pair not needed but initialize for safety
         r1_zbl_pair(:,:) = 0.d0
         r2_zbl_pair(:,:) = 0.d0
         rr_k2b_pair(:,:) = r_cut_pair_in(:,:)

         ! Update global scalar for backward compatibility
         if (r_cut_in < 0.d0) then
            r_cut_in = minval(r_cut_pair_in)
         end if
      end if

      ! =====================================================================
      ! 4. Validate pair-specific parameters
      ! =====================================================================
      if (zbl_potential) then
         call validate_pair_cutoffs(nsp)
      end if

      ! =====================================================================
      ! 5. Log the pair-specific table
      ! =====================================================================
      call log_pair_cutoff_table(nsp)

   end subroutine setup_pair_cutoffs


   !----------------------------------------------------------------------------
   ! apply_pair_override
   !
   ! Parse an override string of the form "El1 El2 val ; El1 El2 val ; ..."
   ! and apply to the pair array. Symmetric: A-B and B-A both set.
   !----------------------------------------------------------------------------
   subroutine apply_pair_override(override_str, pair_array, nsp, label)
      use module_chemical_species, only: fix_ch_elements
      use mld_logger
      use mld_string
      use mld_mpi, only: mld_mpi_abort

      implicit none
      character(len=*), intent(in) :: override_str, label
      integer, intent(in) :: nsp
      real(kind_double), dimension(nsp, nsp), intent(inout) :: pair_array

      character(len=5000) :: work_str
      character(len=100)  :: token
      character(len=2)    :: el1, el2
      real(kind_double)   :: val
      integer :: ipos, isemi, idx1, idx2, ios
      logical :: found1, found2

      work_str = adjustl(override_str)

      do while (len_trim(work_str) > 0)
         ! find semicolon separator
         isemi = index(work_str, ';')
         if (isemi > 0) then
            token = adjustl(work_str(1:isemi-1))
            work_str = adjustl(work_str(isemi+1:))
         else
            token = adjustl(work_str)
            work_str = ''
         end if

         if (len_trim(token) == 0) cycle

         ! parse "El1 El2 val"
         read(token, *, iostat=ios) el1, el2, val
         if (ios /= 0) then
            call mld_mpi_abort('pair_cutoffs: cannot parse "'//label//'" entry: "'//trim(token)//'"'// &
               ' Expected format: "El1 El2 value" (e.g., "Fe W 1.5")')
         end if

         ! find species indices
         found1 = .false.
         found2 = .false.
         do ipos = 1, nsp
            if (fix_ch_elements(ipos) == el1) then
               idx1 = ipos
               found1 = .true.
            end if
            if (fix_ch_elements(ipos) == el2) then
               idx2 = ipos
               found2 = .true.
            end if
         end do

         if (.not. found1) then
            call mld_mpi_abort('pair_cutoffs: element "'//trim(el1)//'" in "'//label// &
               '" is not present in chemical_elements. Check your input.')
         end if
         if (.not. found2) then
            call mld_mpi_abort('pair_cutoffs: element "'//trim(el2)//'" in "'//label// &
               '" is not present in chemical_elements. Check your input.')
         end if

         ! apply symmetrically
         pair_array(idx1, idx2) = val
         pair_array(idx2, idx1) = val

         call log_info('Pair override ('//label//'): '//trim(el1)//' - '//trim(el2)//' = '//vtoa(val)//' Ang')
      end do

   end subroutine apply_pair_override


   !----------------------------------------------------------------------------
   ! validate_pair_cutoffs
   !
   ! Verify the strict ordering r1_zbl_pair < r_cut_pair_in < r2_zbl_pair
   ! for all element pairs when ZBL is active.
   !----------------------------------------------------------------------------
   subroutine validate_pair_cutoffs(nsp)
      use module_neigh_local, only: r_cut_pair_in
      use module_zbl, only: r1_zbl_pair, r2_zbl_pair
      use module_chemical_species, only: fix_ch_elements
      use mld_logger
      use mld_mpi, only: mld_mpi_abort
      use mld_string

      implicit none
      integer, intent(in) :: nsp
      integer :: ii, jj

      do ii = 1, nsp
         do jj = ii, nsp
            if (r1_zbl_pair(ii,jj) >= r_cut_pair_in(ii,jj)) then
               call log_critical('pair_cutoffs: r1_zbl_pair('//trim(fix_ch_elements(ii))//','// &
                  trim(fix_ch_elements(jj))//') = '//vtoa(r1_zbl_pair(ii,jj))// &
                  ' >= r_cut_pair_in = '//vtoa(r_cut_pair_in(ii,jj)))
               call mld_mpi_abort('pair_cutoffs: r1_zbl must be < r_cut_in for all element pairs.')
            end if
            if (r_cut_pair_in(ii,jj) >= r2_zbl_pair(ii,jj)) then
               call log_critical('pair_cutoffs: r_cut_pair_in('//trim(fix_ch_elements(ii))//','// &
                  trim(fix_ch_elements(jj))//') = '//vtoa(r_cut_pair_in(ii,jj))// &
                  ' >= r2_zbl_pair = '//vtoa(r2_zbl_pair(ii,jj)))
               call mld_mpi_abort('pair_cutoffs: r_cut_in must be < r2_zbl for all element pairs.')
            end if
            if (r1_zbl_pair(ii,jj) >= r2_zbl_pair(ii,jj)) then
               call log_critical('pair_cutoffs: r1_zbl_pair('//trim(fix_ch_elements(ii))//','// &
                  trim(fix_ch_elements(jj))//') = '//vtoa(r1_zbl_pair(ii,jj))// &
                  ' >= r2_zbl_pair = '//vtoa(r2_zbl_pair(ii,jj)))
               call mld_mpi_abort('pair_cutoffs: r1_zbl must be < r2_zbl for all element pairs.')
            end if
         end do
      end do

   end subroutine validate_pair_cutoffs


   !----------------------------------------------------------------------------
   ! log_pair_cutoff_table
   !
   ! Print a summary table of pair-specific cutoff parameters.
   !----------------------------------------------------------------------------
   subroutine log_pair_cutoff_table(nsp)
      use module_neigh_local, only: r_cut_pair_in
      use module_zbl, only: zbl_potential, r1_zbl_pair, r2_zbl_pair
      use module_chemical_species, only: fix_ch_elements
      use mld_logger
      use mld_string

      implicit none
      integer, intent(in) :: nsp
      integer :: ii, jj
      character(len=:), allocatable :: msg

      msg = 'Pair-specific cutoff parameters:'//char(10)
      if (zbl_potential) then
         msg = msg//'  Pair        r_cut_in   r1_zbl     r2_zbl'//char(10)
         msg = msg//'  ----------  ---------  ---------  ---------'//char(10)
         do ii = 1, nsp
            do jj = ii, nsp
               msg = msg//'  '//fix_ch_elements(ii)//' - '//fix_ch_elements(jj)// &
                  '     '//vtoa(r_cut_pair_in(ii,jj))// &
                  '  '//vtoa(r1_zbl_pair(ii,jj))// &
                  '  '//vtoa(r2_zbl_pair(ii,jj))//char(10)
            end do
         end do
      else
         msg = msg//'  Pair        r_cut_in'//char(10)
         msg = msg//'  ----------  ---------'//char(10)
         do ii = 1, nsp
            do jj = ii, nsp
               msg = msg//'  '//fix_ch_elements(ii)//' - '//fix_ch_elements(jj)// &
                  '     '//vtoa(r_cut_pair_in(ii,jj))//char(10)
            end do
         end do
      end if

      call log_info(msg)

   end subroutine log_pair_cutoff_table


   !----------------------------------------------------------------------------
   ! write_pair_cutoffs_xml
   !
   ! Write an XML file with all cutoff and ZBL parameters for LAMMPS.
   !----------------------------------------------------------------------------
   subroutine write_pair_cutoffs_xml()
      use module_neigh_local, only: r_cut, r_cut_width, r_cut_width_in, type_fcut, r_cut_pair_in
      use module_zbl, only: zbl_potential, zbl_type, r1_zbl, r2_zbl, &
                            r1_zbl_pair, r2_zbl_pair, rr_k2b_pair, params_k2b_to_zbl
      use module_chemical_species, only: fix_no_of_elements, fix_ch_elements, &
                                         fix_Z_elements, size_species_half
      use ml_in_ndm_module, only: rangml
      use mld_logger
      use mld_string

      implicit none
      integer :: xmlunit, ii, jj, kk, nsp, is
      character(len=200) :: fname

      nsp = fix_no_of_elements

      if (rangml /= 0) return

      fname = 'cutoff_params.xml'
      open(newunit=xmlunit, file=trim(fname), status='replace', action='write')

      write(xmlunit, '(a)') '<?xml version="1.0" encoding="UTF-8"?>'
      write(xmlunit, '(a)') '<cutoff_parameters>'

      ! --- global parameters ---
      write(xmlunit, '(2x,a)') '<global'
      write(xmlunit, '(4x,a,e16.8,a)') 'r_cut="', r_cut, '"'
      write(xmlunit, '(4x,a,e16.8,a)') 'r_cut_width="', r_cut_width, '"'
      write(xmlunit, '(4x,a,i0,a)')    'type_fcut="', type_fcut, '"'
      write(xmlunit, '(4x,a,e16.8,a)') 'r_cut_width_in="', r_cut_width_in, '"'
      if (zbl_potential) then
         write(xmlunit, '(4x,a,i0,a)')    'zbl_type="', zbl_type, '"'
      end if
      write(xmlunit, '(2x,a)') '/>'

      ! --- pair cutoffs ---
      write(xmlunit, '(2x,a,i0,a)', advance='no') '<pair_cutoffs n_species="', nsp, '" elements="'
      do ii = 1, nsp
         if (ii > 1) write(xmlunit, '(a)', advance='no') ' '
         write(xmlunit, '(a)', advance='no') trim(fix_ch_elements(ii))
      end do
      write(xmlunit, '(a)') '">'

      do ii = 1, nsp
         do jj = ii, nsp
            write(xmlunit, '(4x,a)', advance='no') '<pair'
            write(xmlunit, '(a,a,a)', advance='no') ' el1="', trim(fix_ch_elements(ii)), '"'
            write(xmlunit, '(a,a,a)', advance='no') ' el2="', trim(fix_ch_elements(jj)), '"'
            write(xmlunit, '(a,e16.8,a)', advance='no') ' r_cut_in="', r_cut_pair_in(ii,jj), '"'
            if (zbl_potential) then
               write(xmlunit, '(a,e16.8,a)', advance='no') ' r1_zbl="', r1_zbl_pair(ii,jj), '"'
               write(xmlunit, '(a,e16.8,a)', advance='no') ' r2_zbl="', r2_zbl_pair(ii,jj), '"'
               write(xmlunit, '(a,e16.8,a)', advance='no') ' rr_k2b="', rr_k2b_pair(ii,jj), '"'
            end if
            write(xmlunit, '(a)') ' />'
         end do
      end do
      write(xmlunit, '(2x,a)') '</pair_cutoffs>'

      ! --- ZBL bridge parameters (if present) ---
      if (zbl_potential .and. allocated(params_k2b_to_zbl)) then
         write(xmlunit, '(2x,a,i0,a)') '<zbl_bridge zbl_type="', zbl_type, '">'
         do is = 1, size(params_k2b_to_zbl, 2)
            write(xmlunit, '(4x,a,i0,a)', advance='no') '<bridge_params species_pair="', is, '"'
            do kk = 1, size(params_k2b_to_zbl, 1)
               write(xmlunit, '(a,i0,a,e20.12,a)', advance='no') ' p', kk, '="', params_k2b_to_zbl(kk, is), '"'
            end do
            write(xmlunit, '(a)') ' />'
         end do
         write(xmlunit, '(2x,a)') '</zbl_bridge>'
      end if

      ! --- ZBL constants (for reference) ---
      if (zbl_potential) then
         write(xmlunit, '(2x,a)') '<zbl_constants>'
         do ii = 1, nsp
            write(xmlunit, '(4x,a,a,a,e16.8,a)') &
               '<element symbol="', trim(fix_ch_elements(ii)), '" Z="', fix_Z_elements(ii), '" />'
         end do
         write(xmlunit, '(2x,a)') '</zbl_constants>'
      end if

      write(xmlunit, '(a)') '</cutoff_parameters>'
      close(xmlunit)

      call log_info('Pair cutoff parameters written to '//trim(fname))

   end subroutine write_pair_cutoffs_xml

end module module_pair_cutoffs
