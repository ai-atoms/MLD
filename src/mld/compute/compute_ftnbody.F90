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

#include "../../MLD_MACROS.INC"

!===============================================================================
! module_ftnbody_potio
!
! Persistence of the ftnbody random-feature descriptor data (random.xml).
! ftnbody is a random-sampling descriptor: its frequencies (omega) and phases
! depend on the training data (via find_best_length), so they cannot be
! regenerated identically at prediction/MD time. They are therefore written to
! random.xml during the fit and read back when a saved potential is loaded.
! (Low-rank embeddings for chem_mode=2 are a planned addition to this file.)
!===============================================================================
module module_ftnbody_potio
  use module_kind_variables, only: kind_double
  implicit none
  private
  public :: write_ftnbody_random_xml, read_ftnbody_random_xml

contains

  subroutine write_ftnbody_random_xml(fname)
    use module_ftnbody, only: ftnbody_omega_2b, ftnbody_omega_3b, ftnbody_omega_4b, ftnbody_omega_5b, &
                              ftnbody_phase_random_2b, ftnbody_phase_random_3b, ftnbody_phase_random_4b, &
                              ftnbody_phase_random_5b, ftnbody_chem_mode, &
                              ftnbody_chem_rank_n, ftnbody_hash_channels_n, &
                              ftnbody_chem_center, ftnbody_chem_neigh, ftnbody_model, &
                              spip_degree_n
    use module_body_desc, only: l_body_order
    use module_chemical_species, only: fix_no_of_elements
    character(len=*), intent(in) :: fname
    integer :: u, u_i
    open(newunit=u, file=trim(fname), status='replace', action='write')
    write(u,'(a)') '<?xml version="1.0"?>'
    write(u,'(a)') '<ftnbody_random_features>'
    write(u,'(2x,a,i0,a,i0,a,a,a)') '<meta chem_mode="', ftnbody_chem_mode, &
                                '" n_species="', fix_no_of_elements, &
                                '" model="', trim(ftnbody_model), '"/>'
    ! chemistry layout (skipped by the numeric reader; parsed by the LAMMPS C++ reader)
    write(u,'(2x,a)', advance='no') '<chem chem_rank="'
    do u_i = 2, 5; write(u,'(i0,1x)', advance='no') ftnbody_chem_rank_n(u_i); end do
    write(u,'(a)', advance='no') '" hash_channels="'
    do u_i = 2, 5; write(u,'(i0,1x)', advance='no') ftnbody_hash_channels_n(u_i); end do
    ! spip: the orbit tables are regenerated deterministically from the degree,
    ! so only spip_degree needs to travel with the potential
    write(u,'(a)', advance='no') '" spip_degree="'
    do u_i = 2, 5; write(u,'(i0,1x)', advance='no') spip_degree_n(u_i); end do
    write(u,'(a)') '"/>'
    if (l_body_order(2)) call write_order_block(u, 2, ftnbody_omega_2b, ftnbody_phase_random_2b)
    if (l_body_order(3)) call write_order_block(u, 3, ftnbody_omega_3b, ftnbody_phase_random_3b)
    if (l_body_order(4)) call write_order_block(u, 4, ftnbody_omega_4b, ftnbody_phase_random_4b)
    if (l_body_order(5)) call write_order_block(u, 5, ftnbody_omega_5b, ftnbody_phase_random_5b)
    ! mode 2 (low-rank): the role embeddings are data-/RNG-dependent and must be
    ! persisted so LAMMPS (no Fortran RNG) can reconstruct the chemistry weights.
    if (ftnbody_chem_mode == 2 .and. allocated(ftnbody_chem_center)) then
      write(u,'(2x,a,i0,a,i0,a,i0,a,i0,a)') '<embeddings n_species="', size(ftnbody_chem_center,1), &
            '" nchan="', size(ftnbody_chem_center,2), '" nroles="', size(ftnbody_chem_neigh,3), &
            '" max_order="', size(ftnbody_chem_center,3), '">'
      call write_real_dump(u, '<chem_center>', '</chem_center>', ftnbody_chem_center, size(ftnbody_chem_center))
      call write_real_dump(u, '<chem_neigh>',  '</chem_neigh>',  ftnbody_chem_neigh,  size(ftnbody_chem_neigh))
      write(u,'(2x,a)') '</embeddings>'
    end if
    write(u,'(a)') '</ftnbody_random_features>'
    close(u)
  end subroutine write_ftnbody_random_xml

  ! dump a real array (column-major sequence order), one value per line, wrapped in tags
  subroutine write_real_dump(u, tag_open, tag_close, arr, n)
    integer, intent(in) :: u, n
    character(len=*), intent(in) :: tag_open, tag_close
    real(kind_double), intent(in) :: arr(*)
    integer :: k
    write(u,'(4x,a)') tag_open
    do k = 1, n
      write(u,'(6x,ES25.15E3)') arr(k)
    end do
    write(u,'(4x,a)') tag_close
  end subroutine write_real_dump

  subroutine write_order_block(u, n, omega, phase)
    integer, intent(in) :: u, n
    real(kind_double), intent(in) :: omega(:,:), phase(:)
    integer :: i, j
    write(u,'(2x,a,i0,a,i0,a,i0,a)') '<order n="', n, '" rows="', size(omega,1), &
                                     '" cols="', size(omega,2), '">'
    write(u,'(4x,a)') '<omega>'
    do j = 1, size(omega,2)
      do i = 1, size(omega,1)
        write(u,'(6x,ES25.15E3)') omega(i,j)
      end do
    end do
    write(u,'(4x,a)') '</omega>'
    write(u,'(4x,a)') '<phase>'
    do i = 1, size(phase,1)
      write(u,'(6x,ES25.15E3)') phase(i)
    end do
    write(u,'(4x,a)') '</phase>'
    write(u,'(2x,a)') '</order>'
  end subroutine write_order_block

  subroutine read_ftnbody_random_xml(fname)
    use module_ftnbody, only: ftnbody_omega_2b, ftnbody_omega_3b, ftnbody_omega_4b, ftnbody_omega_5b, &
                              ftnbody_phase_random_2b, ftnbody_phase_random_3b, ftnbody_phase_random_4b, &
                              ftnbody_phase_random_5b, ftnbody_chem_mode, &
                              ftnbody_chem_center, ftnbody_chem_neigh
    use module_body_desc, only: l_body_order
    use mld_logger
    character(len=*), intent(in) :: fname
    integer :: u
    logical :: ok
    inquire(file=trim(fname), exist=ok)
    if (.not. ok) then
      call log_critical("read_ftnbody_random_xml: missing "//trim(fname)// &
                        " (needed to reproduce the ftnbody descriptor for MD/prediction)")
      return
    end if
    open(newunit=u, file=trim(fname), status='old', action='read')
    ! same active-order order as the writer
    if (l_body_order(2)) call read_order_block(u, ftnbody_omega_2b, ftnbody_phase_random_2b)
    if (l_body_order(3)) call read_order_block(u, ftnbody_omega_3b, ftnbody_phase_random_3b)
    if (l_body_order(4)) call read_order_block(u, ftnbody_omega_4b, ftnbody_phase_random_4b)
    if (l_body_order(5)) call read_order_block(u, ftnbody_omega_5b, ftnbody_phase_random_5b)
    ! mode 2: load the persisted low-rank role embeddings (same column-major order
    ! as the writer); the arrays were allocated in init_ftnbody.
    if (ftnbody_chem_mode == 2 .and. allocated(ftnbody_chem_center)) then
      call read_n_reals(u, size(ftnbody_chem_center), ftnbody_chem_center)
      call read_n_reals(u, size(ftnbody_chem_neigh),  ftnbody_chem_neigh)
    end if
    close(u)
    call log_info("read_ftnbody_random_xml: ftnbody random features loaded from "//trim(fname))
  end subroutine read_ftnbody_random_xml

  subroutine read_order_block(u, omega, phase)
    integer, intent(in) :: u
    real(kind_double), intent(inout) :: omega(:,:), phase(:)
    ! omega is filled column-major (matches the writer's j-outer/i-inner order)
    call read_n_reals(u, size(omega), omega)
    call read_n_reals(u, size(phase), phase)
  end subroutine read_order_block

  ! read n reals, one per line, skipping any XML tag or blank line
  subroutine read_n_reals(u, n, arr)
    integer, intent(in) :: u, n
    real(kind_double), intent(out) :: arr(*)
    integer :: k, ios
    character(len=512) :: line
    k = 0
    do while (k < n)
      read(u,'(a)', iostat=ios) line
      if (ios /= 0) exit
      if (index(line,'<') > 0) cycle
      if (len_trim(line) == 0) cycle
      k = k + 1
      read(line,*) arr(k)
    end do
  end subroutine read_n_reals

end module module_ftnbody_potio


  subroutine init_ftnbody
    use module_kind_variables, only: kind_double
    use module_body_desc, only: l_body_order, dim_desc_body 
    use module_ftnbody, only: ftnbody_dim, dim_qbody, dim_rff, length_rff, ftnbody_phase_random_2b, &
            ftnbody_phase_random_3b, ftnbody_phase_random_4b, ftnbody_phase_random_5b, &
            ftnbody_omega_2b, ftnbody_omega_3b, ftnbody_omega_4b, transpose_ftnbody_omega_4b, ftnbody_omega_5b, &
            covar_matrix_2b, covar_matrix_3b, covar_matrix_4b, covar_matrix_5b, &
            length_rff_2b, length_rff_3b, length_rff_4b, length_rff_5b , count_2b, count_3b, count_4b, count_5b, &
            find_best_length, mean_2b, mean_3b, mean_4b, mean_5b, &
            ftnbody_chem_mode, ftnbody_chem_rank_n, ftnbody_hash_channels_n, ftnbody_n_channels, &
            ftnbody_chem_center, ftnbody_chem_neigh, MAX_ORDER_FTNBODY, ftnbody_read_from_potential, &
            ftnbody_model_id, FTNBODY_MODEL_GRAMM, FTNBODY_MODEL_CPIP, FTNBODY_MODEL_SPIP, &
            ftnbody_dimq, spip_tab
    ! use module_kernel, only: np_kernel_full
    use module_kernel, only:  krff_type
    use module_chemical_species, only: fix_no_of_elements
    use module_ftnbody_potio, only: read_ftnbody_random_xml
    use mld_mpi, only: mld_critical_abort
    use mld_logger
    use module_sample_rand_ker, only: get_rff_val_sigma
    implicit none
    real(kind_double)    :: mu
    integer :: ii, icnt
    integer :: nord, nchan, irk, isp, nseed
    integer, allocatable :: iseed(:)
      
    _NAMECURRENT_("init_ftnbody")
    _MLD_BEGIN_

    mu = 0.d0

    dim_desc_body(:)=0

    ! geometry model (docs/MLT5 sec. 2.2.4 and 2.2.6).  The per-order RFF input
    ! dimension dim_qbody(n) = ftnbody_dimq(n) depends on the model:
    !   poly/gramm : d_n = n(n-1)/2 compact / ordered Gram coordinates
    !   cpip       : CPIP_MG_N(n) geometric (+ CPIP_MC_N(n) colored, chem mode 2)
    !   spip       : M_n monomial orbits up to spip_degree_nbody (built here)
    ! All allocations below are expressed in dim_qbody(n), so they are common.
    select case (ftnbody_model_id)
    case (FTNBODY_MODEL_GRAMM)
      call init_ftnbody_gram_perms
      call log_info("ML: FT-nBody geometry model gramm: ordered Gram coordinates "// &
                    "(u_1..u_m, c_pq), exact S_{n-1} orbit averaging of the Fourier features")
    case (FTNBODY_MODEL_CPIP)
      if (ftnbody_chem_mode /= 0 .and. ftnbody_chem_mode /= 2) &
        call mld_critical_abort("init_ftnbody: ftnbody_model='cpip' supports only chem modes 0 and 2")
      call log_info("ML: FT-nBody geometry model cpip: compact mixed invariant set "// &
                    "(radial/angular/mixed moments; chem mode 2 adds colored moments to the RFF input)")
    case (FTNBODY_MODEL_SPIP)
      call init_ftnbody_gram_perms   ! S_{n-1} permutation tables (shared with gramm)
      call init_ftnbody_spip_orbits  ! monomial orbit tables up to spip_degree_nbody
      call log_info("ML: FT-nBody geometry model spip: systematic S_{n-1} monomial-orbit "// &
                    "coordinates up to spip_degree_nbody")
    case default
      call log_info("ML: FT-nBody geometry model poly: compact invariant coordinates")
    end select

    ! find_best_length = .true.
    if (.not.(find_best_length)) then
      call log_info("ML: FT-nBody warning: without finding best length!!!")
    end if

    if (l_body_order(2)) then
      ! set the F 
      ! dim_rff(2) =  np_kernel_full
      ! size of the q vectors (model dependent, see ftnbody_dimq) ...
      dim_qbody(2) = ftnbody_dimq(2)
      ! true only for one typer of atoms  ... 
      dim_desc_body(2) = dim_rff(2)  
      !TODOspec for multispecies: 

      count_2b = 0.d0
      if (allocated(covar_matrix_2b)) deallocate (covar_matrix_2b); allocate (covar_matrix_2b(dim_qbody(2), dim_qbody(2)))
      if (allocated(length_rff_2b)) deallocate (length_rff_2b); allocate (length_rff_2b(dim_qbody(2)))
      if (allocated(mean_2b)) deallocate (mean_2b); allocate (mean_2b(dim_qbody(2)))
      covar_matrix_2b(:,:) = 0.d0
      mean_2b(:)=0.d0 

      if (allocated(ftnbody_omega_2b)) deallocate (ftnbody_omega_2b); allocate (ftnbody_omega_2b(dim_qbody(2), dim_desc_body(2)))
      if (allocated(ftnbody_phase_random_2b)) deallocate (ftnbody_phase_random_2b); allocate (ftnbody_phase_random_2b(dim_desc_body(2)))
      if (.not.(find_best_length)) then
        call  get_rff_val_sigma(krff_type, mu, dim_qbody(2), dim_desc_body(2), ftnbody_omega_2b, ftnbody_phase_random_2b)
        ftnbody_omega_2b = ftnbody_omega_2b / (length_rff(2)*sqrt(2.d0))
      end if 

    end if   


    if (l_body_order(3)) then
      ! set the F 
      ! dim_rff(3) =  np_kernel_full
      ! size of the q vectors (model dependent, see ftnbody_dimq) ...
      dim_qbody(3) = ftnbody_dimq(3)
      ! true only for one typer of atoms  ... 
      dim_desc_body(3) = dim_rff(3)  
      !TODOspec for multispecies: 

      count_3b = 0.d0
      if (allocated(covar_matrix_3b)) deallocate (covar_matrix_3b); allocate (covar_matrix_3b(dim_qbody(3), dim_qbody(3)))
      if (allocated(length_rff_3b)) deallocate (length_rff_3b); allocate (length_rff_3b(dim_qbody(3)))
      if (allocated(mean_3b)) deallocate (mean_3b); allocate (mean_3b(dim_qbody(3)))
      covar_matrix_3b(:,:) = 0.d0
      mean_3b(:)=0.d0 

      if (allocated(ftnbody_omega_3b)) deallocate (ftnbody_omega_3b); allocate (ftnbody_omega_3b(dim_qbody(3), dim_desc_body(3)))
      if (allocated(ftnbody_phase_random_3b)) deallocate (ftnbody_phase_random_3b); allocate (ftnbody_phase_random_3b(dim_desc_body(3)))

      if (.not.(find_best_length)) then
        call  get_rff_val_sigma(krff_type, mu, dim_qbody(3), dim_desc_body(3), ftnbody_omega_3b, ftnbody_phase_random_3b)
        ftnbody_omega_3b = ftnbody_omega_3b / (length_rff(3)*sqrt(2.d0))
      end if

    end if   

    if (l_body_order(4)) then
      ! set the F 
      ! dim_rff(4) =  np_kernel_full
      ! size of the q vectors (model dependent, see ftnbody_dimq) ...
      dim_qbody(4) = ftnbody_dimq(4)
      ! true only for one typer of atoms  ... 
      dim_desc_body(4) = dim_rff(4)  
      !TODOspec for multispecies: 

      count_4b = 0.d0
      if (allocated(covar_matrix_4b)) deallocate (covar_matrix_4b); allocate (covar_matrix_4b(dim_qbody(4), dim_qbody(4)))
      if (allocated(length_rff_4b)) deallocate (length_rff_4b); allocate (length_rff_4b(dim_qbody(4)))
      if (allocated(mean_4b)) deallocate (mean_4b); allocate (mean_4b(dim_qbody(4)))
      covar_matrix_4b(:,:) = 0.d0
      mean_4b(:)=0.d0 

      if (allocated(ftnbody_omega_4b)) deallocate (ftnbody_omega_4b); allocate (ftnbody_omega_4b(dim_qbody(4), dim_desc_body(4)))
      if (allocated(transpose_ftnbody_omega_4b)) deallocate (transpose_ftnbody_omega_4b); allocate (transpose_ftnbody_omega_4b(dim_desc_body(4),  dim_qbody(4)))
      if (allocated(ftnbody_phase_random_4b)) deallocate (ftnbody_phase_random_4b); allocate (ftnbody_phase_random_4b(dim_desc_body(4)))

      if (.not.(find_best_length)) then
        call  get_rff_val_sigma(krff_type, mu, dim_qbody(4), dim_desc_body(4), ftnbody_omega_4b, ftnbody_phase_random_4b)
        ftnbody_omega_4b = ftnbody_omega_4b / (length_rff(4)*sqrt(2.d0))
      end if

      transpose_ftnbody_omega_4b = transpose(ftnbody_omega_4b)
        
    end if   


    if (l_body_order(5)) then
      ! set the F 
      ! dim_rff(5) =  np_kernel_full
      ! size of the q vectors (model dependent, see ftnbody_dimq) ...
      dim_qbody(5) = ftnbody_dimq(5)
      ! true only for one typer of atoms  ... 
      dim_desc_body(5) = dim_rff(5)  
      !TODOspec for multispecies: 

      count_5b = 0.d0
      if (allocated(covar_matrix_5b)) deallocate (covar_matrix_5b); allocate (covar_matrix_5b(dim_qbody(5), dim_qbody(5)))
      if (allocated(length_rff_5b)) deallocate (length_rff_5b); allocate (length_rff_5b(dim_qbody(5)))
      if (allocated(mean_5b)) deallocate (mean_5b); allocate (mean_5b(dim_qbody(5)))
      covar_matrix_5b(:,:) = 0.d0
      mean_5b(:)=0.d0 

      if (allocated(ftnbody_omega_5b)) deallocate (ftnbody_omega_5b); allocate (ftnbody_omega_5b(dim_qbody(5), dim_desc_body(5)))
      if (allocated(ftnbody_phase_random_5b)) deallocate (ftnbody_phase_random_5b); allocate (ftnbody_phase_random_5b(dim_desc_body(5)))

      if (.not.(find_best_length)) then
        call  get_rff_val_sigma(krff_type, mu, dim_qbody(5), dim_desc_body(5), ftnbody_omega_5b, ftnbody_phase_random_5b)
        ftnbody_omega_5b = ftnbody_omega_5b / (length_rff(5)*sqrt(2.d0))
      end if

    end if   

    ! ---------------------------------------------------------------------
    ! Multispecies chemistry modes (see docs/perspective_ftnbody.md).
    ! Geometry (omega/phase) stays sized by dim_rff(n); chemistry expands the
    ! descriptor block: dim_desc_body(n) = dim_rff(n) * ftnbody_n_channels(n).
    !   mode 0 -> 1 channel (single species, unchanged)
    !   mode 1 -> N_chem(n,S) exact channels
    !   mode 2 -> ftnbody_chem_rank_n(n) low-rank channels
    !   mode 3 -> ftnbody_hash_channels_n(n) hashed channels
    ! ---------------------------------------------------------------------
    ftnbody_n_channels(:) = 1
    if (ftnbody_chem_mode /= 0) then
      do nord = 2, MAX_ORDER_FTNBODY
        if (.not. l_body_order(nord)) cycle
        if (nord > 5) then
          call log_critical("init_ftnbody: ftnbody_chem_mode>0 is implemented only up to 5-body; order " &
                            //vtoa(nord)//" with a chemistry mode is not available")
          call mld_critical_abort("init_ftnbody: chem mode for order>5 not implemented")
        end if
        select case (ftnbody_chem_mode)
        case (1)   ! exact channels: N_chem(n,S) = S * C(S+n-2, n-1)
          select case (nord)
          case (2); nchan = fix_no_of_elements * fix_no_of_elements                              ! S*S
          case (3); nchan = fix_no_of_elements * fix_no_of_elements * (fix_no_of_elements+1) / 2 ! S*S*(S+1)/2
          case (4); nchan = fix_no_of_elements * fix_no_of_elements * (fix_no_of_elements+1) &
                                                * (fix_no_of_elements+2) / 6                      ! S^2*(S+1)*(S+2)/6
          case (5); nchan = fix_no_of_elements * fix_no_of_elements * (fix_no_of_elements+1) &
                                                * (fix_no_of_elements+2) * (fix_no_of_elements+3) / 24  ! S^2*(S+1)*(S+2)*(S+3)/24
          end select
        case (2)   ! low-rank embedding
          if (ftnbody_chem_rank_n(nord) <= 0) &
            call mld_critical_abort("init_ftnbody: ftnbody_chem_rank must be > 0 for chem mode 2 at this order")
          nchan = ftnbody_chem_rank_n(nord)
        case (3)   ! hashed channels
          if (ftnbody_hash_channels_n(nord) <= 0) &
            call mld_critical_abort("init_ftnbody: ftnbody_hash_channels must be > 0 for chem mode 3 at this order")
          nchan = ftnbody_hash_channels_n(nord)
        case default
          call mld_critical_abort("init_ftnbody: unknown ftnbody_chem_mode (must be 0,1,2,3)")
        end select
        ftnbody_n_channels(nord) = nchan
        dim_desc_body(nord) = dim_rff(nord) * nchan
        call log_info("ML: FT-nBody chem mode "//vtoa(ftnbody_chem_mode)//" order "//vtoa(nord) &
                     //" -> "//vtoa(nchan)//" channels, dim_desc_body "//vtoa(dim_desc_body(nord)))
      end do

      ! Low-rank role embeddings (mode 2): fixed deterministic-random, column-normalized.
      ! NOTE: these embeddings must be written to / read from the potential file for
      !       prediction (see todo: potential-file export). Same-run train+test is consistent.
      if (ftnbody_chem_mode == 2) then
        nchan = maxval(ftnbody_chem_rank_n(:))
        if (allocated(ftnbody_chem_center)) deallocate(ftnbody_chem_center)
        if (allocated(ftnbody_chem_neigh))  deallocate(ftnbody_chem_neigh)
        allocate(ftnbody_chem_center(fix_no_of_elements, nchan, MAX_ORDER_FTNBODY))
        allocate(ftnbody_chem_neigh (fix_no_of_elements, nchan, MAX_ORDER_FTNBODY-1, MAX_ORDER_FTNBODY))
        ftnbody_chem_center(:,:,:)   = 0.d0
        ftnbody_chem_neigh (:,:,:,:) = 0.d0
        ! At prediction/MD/LAMMPS time the embeddings come from random.xml (read
        ! below); only generate them during a fresh fit.
        if (.not. ftnbody_read_from_potential) then
          call random_seed(size=nseed); allocate(iseed(nseed))
          do ii = 1, nseed
            iseed(ii) = 1234567 + 89*ii
          end do
          call random_seed(put=iseed)
          do nord = 2, 5
            if (.not. l_body_order(nord)) cycle
            do irk = 1, ftnbody_chem_rank_n(nord)
              do isp = 1, fix_no_of_elements
                call random_number(mu); ftnbody_chem_center(isp, irk, nord) = mu - 0.5d0
                do ii = 1, nord - 1   ! one role embedding per neighbour slot
                  call random_number(mu); ftnbody_chem_neigh (isp, irk, ii, nord) = mu - 0.5d0
                end do
              end do
            end do
          end do
          deallocate(iseed)
        end if
      end if
    end if

    icnt = 0
    do ii = 1, size(l_body_order,1)
      if (l_body_order(ii)) icnt = icnt + dim_desc_body(ii)
    end do

    ftnbody_dim = icnt

    do ii = 1, size(dim_desc_body,1)
      if (l_body_order(ii)) then 
      call log_info("ML: FT-nBody with "//vtoa(ii)//"-body has the size "//vtoa(dim_desc_body(ii)))
      end if 
    end do   
    call log_info("ML: FT-nBody with total size "// vtoa(ftnbody_dim))

    ! Prediction / MD from a saved potential: the random features cannot be
    ! re-estimated identically, so load them from random.xml (init_length is
    ! skipped). The descriptor arrays were just allocated above.
    if (ftnbody_read_from_potential) then
      call read_ftnbody_random_xml('random.xml')
    end if

    _MLD_END_
  end subroutine init_ftnbody


  !-----------------------------------------------------------------------------
  ! Gram model (ftnbody_model='gramm'): build, once, the neighbour-permutation
  ! tables of every body order.  For order n with m=n-1 neighbours:
  !   gram_perm(p,ip,n)  = pi(p)  for the ip-th permutation pi of (1..m)
  !   gram_qmap(i,ip,n)  = source index so that (Pi_pi q)_i = q(gram_qmap(i,ip,n))
  ! with the Gram layout q = (u_1..u_m, c_12, c_13, ..., c_{m-1,m}):
  !   (Pi_pi q)_p        = u_{pi(p)}                         (radial slots)
  !   (Pi_pi q)_{pos(p,q)} = c_{pi(p)pi(q)} = q(pos(sorted))  (angular slots)
  !-----------------------------------------------------------------------------
  subroutine init_ftnbody_gram_perms
    use module_ftnbody, only: MAX_ORDER_FTNBODY, gram_nperm, gram_perm, gram_qmap, gram_pair_pos
    implicit none
    integer :: nord, m, ip, i1, i2, i3, i4
    integer :: pv(4)
    integer, parameter, dimension(0:4) :: factorial_m = (/ 1, 1, 2, 6, 24 /)

    gram_nperm(:)     = 1
    gram_perm(:,:,:)  = 0
    gram_qmap(:,:,:)  = 0

    do nord = 2, MAX_ORDER_FTNBODY
      m = nord - 1
      gram_nperm(nord) = factorial_m(m)
      ip = 0
      do i1 = 1, m
        pv(1) = i1
        if (m == 1) then
          ip = ip + 1
          call gram_store_perm(nord, m, ip, pv)
          cycle
        end if
        do i2 = 1, m
          if (i2 == i1) cycle
          pv(2) = i2
          if (m == 2) then
            ip = ip + 1
            call gram_store_perm(nord, m, ip, pv)
            cycle
          end if
          do i3 = 1, m
            if (i3 == i1 .or. i3 == i2) cycle
            pv(3) = i3
            if (m == 3) then
              ip = ip + 1
              call gram_store_perm(nord, m, ip, pv)
              cycle
            end if
            do i4 = 1, m
              if (i4 == i1 .or. i4 == i2 .or. i4 == i3) cycle
              pv(4) = i4
              ip = ip + 1
              call gram_store_perm(nord, m, ip, pv)
            end do
          end do
        end do
      end do
    end do

  contains

    subroutine gram_store_perm(l_nord, l_m, l_ip, l_pv)
      integer, intent(in) :: l_nord, l_m, l_ip, l_pv(4)
      integer :: p, q, a, b
      gram_perm(1:l_m, l_ip, l_nord) = l_pv(1:l_m)
      do p = 1, l_m
        gram_qmap(p, l_ip, l_nord) = l_pv(p)
      end do
      do p = 1, l_m - 1
        do q = p + 1, l_m
          a = min(l_pv(p), l_pv(q))
          b = max(l_pv(p), l_pv(q))
          gram_qmap(gram_pair_pos(l_m, p, q), l_ip, l_nord) = gram_pair_pos(l_m, a, b)
        end do
      end do
    end subroutine gram_store_perm

  end subroutine init_ftnbody_gram_perms


  !-----------------------------------------------------------------------------
  ! sPIP model (ftnbody_model='spip', docs/MLT5 sec. "Level C: sPIP"): build,
  ! once, the monomial orbit tables of every active body order.  For order n
  ! with m = n-1 neighbours the primitive coordinates are
  !   x = (u_1..u_m, c_12, c_13, ..., c_{m-1,m}),  dimx = m + C(m,2),
  ! and every monomial x^e with 1 <= |e| <= spip_degree_n(n) belongs to one
  ! S_m orbit.  One representative exponent vector is kept per orbit together
  ! with the collapsed distinct-member list, so that at run time
  !   P_alpha(x) = sum_k mono_coef(k) * x^mono_expo(:,k)     (Reynolds average)
  ! without any group operation inside the atomic loop.  Requires the
  ! permutation tables of init_ftnbody_gram_perms (gram_qmap).
  !-----------------------------------------------------------------------------
  subroutine init_ftnbody_spip_orbits
    use module_kind_variables, only: kind_double
    use module_body_desc, only: l_body_order
    use module_ftnbody, only: MAX_ORDER_FTNBODY, spip_tab, spip_degree_n, &
                              gram_nperm, gram_qmap
    use mld_mpi, only: mld_critical_abort
    use mld_logger
    implicit none
    integer, parameter :: MAXX = 10, MAXPER = 24
    integer :: nord, m, dimx, ndeg, nper, nvec
    integer :: i, ip, iorb, imono, im, nmem, idg
    integer :: e(MAXX), eimg(MAXX)
    integer :: mem_expo(MAXX, MAXPER), mem_mult(MAXPER)
    integer, allocatable :: tmp_rep(:, :), tmp_expo(:, :), tmp_ptr(:)
    real(kind_double), allocatable :: tmp_coef(:)
    logical :: canon, found, done

    do nord = 2, MAX_ORDER_FTNBODY
      if (.not. l_body_order(nord)) cycle
      m    = nord - 1
      dimx = m + m*(m - 1)/2
      ndeg = spip_degree_n(nord)
      nper = gram_nperm(nord)
      if (ndeg < 1) &
        call mld_critical_abort("init_ftnbody_spip_orbits: spip_degree_nbody must be >= 1 for active order " &
                                //vtoa(nord))
      if (ndeg > 10) &
        call mld_critical_abort("init_ftnbody_spip_orbits: spip_degree_nbody > 10 not supported (order " &
                                //vtoa(nord)//")")

      ! number of monomials with 1 <= |e| <= ndeg: C(dimx+ndeg, ndeg) - 1
      nvec = 1
      do idg = 1, ndeg
        nvec = nvec*(dimx + idg)/idg
      end do
      nvec = nvec - 1
      allocate(tmp_rep(dimx, nvec), tmp_expo(dimx, nvec), tmp_ptr(nvec + 1), tmp_coef(nvec))

      iorb  = 0
      imono = 0
      e(1:dimx) = 0
      done = .false.
      odometer: do
        ! advance the exponent vector to the next |e| <= ndeg (mixed radix)
        i = 1
        advance: do
          e(i) = e(i) + 1
          if (sum(e(1:dimx)) <= ndeg) exit advance
          e(i) = 0
          i = i + 1
          if (i > dimx) then
            done = .true.
            exit advance
          end if
        end do advance
        if (done) exit odometer

        ! e is an orbit representative iff no permuted image is lexicographically
        ! greater; the image under pi carries e(i) to slot gram_qmap(i,pi,nord)
        canon = .true.
        do ip = 1, nper
          do i = 1, dimx
            eimg(gram_qmap(i, ip, nord)) = e(i)
          end do
          if (lex_gt(eimg, e, dimx)) then
            canon = .false.
            exit
          end if
        end do
        if (.not. canon) cycle odometer

        ! new orbit: collapse the m! images into distinct members + multiplicity
        nmem = 0
        do ip = 1, nper
          do i = 1, dimx
            eimg(gram_qmap(i, ip, nord)) = e(i)
          end do
          found = .false.
          do im = 1, nmem
            if (all(mem_expo(1:dimx, im) == eimg(1:dimx))) then
              mem_mult(im) = mem_mult(im) + 1
              found = .true.
              exit
            end if
          end do
          if (.not. found) then
            nmem = nmem + 1
            mem_expo(1:dimx, nmem) = eimg(1:dimx)
            mem_mult(nmem) = 1
          end if
        end do

        iorb = iorb + 1
        tmp_rep(1:dimx, iorb) = e(1:dimx)
        tmp_ptr(iorb) = imono + 1
        do im = 1, nmem
          imono = imono + 1
          tmp_expo(1:dimx, imono) = mem_expo(1:dimx, im)
          tmp_coef(imono) = dble(mem_mult(im))/dble(nper)
        end do
      end do odometer
      tmp_ptr(iorb + 1) = imono + 1

      spip_tab(nord)%norb  = iorb
      spip_tab(nord)%nmono = imono
      if (allocated(spip_tab(nord)%orb_ptr))   deallocate(spip_tab(nord)%orb_ptr)
      if (allocated(spip_tab(nord)%mono_expo)) deallocate(spip_tab(nord)%mono_expo)
      if (allocated(spip_tab(nord)%mono_coef)) deallocate(spip_tab(nord)%mono_coef)
      if (allocated(spip_tab(nord)%rep_expo))  deallocate(spip_tab(nord)%rep_expo)
      allocate(spip_tab(nord)%orb_ptr(iorb + 1))
      allocate(spip_tab(nord)%mono_expo(dimx, imono))
      allocate(spip_tab(nord)%mono_coef(imono))
      allocate(spip_tab(nord)%rep_expo(dimx, iorb))
      spip_tab(nord)%orb_ptr(1:iorb + 1)     = tmp_ptr(1:iorb + 1)
      spip_tab(nord)%mono_expo(:, 1:imono)   = tmp_expo(:, 1:imono)
      spip_tab(nord)%mono_coef(1:imono)      = tmp_coef(1:imono)
      spip_tab(nord)%rep_expo(:, 1:iorb)     = tmp_rep(:, 1:iorb)
      deallocate(tmp_rep, tmp_expo, tmp_ptr, tmp_coef)

      call log_info("ML: FT-nBody spip order "//vtoa(nord)//": degree "//vtoa(ndeg) &
                   //" -> "//vtoa(iorb)//" orbit coordinates ("//vtoa(imono)//" member monomials)")
    end do

  contains

    ! lexicographic a > b for exponent vectors of length n
    pure logical function lex_gt(a, b, n) result(gt)
      integer, intent(in) :: a(*), b(*), n
      integer :: k
      gt = .false.
      do k = 1, n
        if (a(k) > b(k)) then
          gt = .true.
          return
        else if (a(k) < b(k)) then
          return
        end if
      end do
    end function lex_gt

  end subroutine init_ftnbody_spip_orbits


  subroutine init_length_ftnbody
    use module_kind_variables, only: kind_double
    use module_body_desc, only: l_body_order
    use module_body_desc, only: dim_desc_body
    use module_ftnbody, only: dim_qbody, dim_rff, init_mode_ftnbody, covar_matrix_2b, covar_matrix_3b, covar_matrix_4b, covar_matrix_5b, count_2b, count_3b, count_4b, count_5b, &
                              length_rff_2b, length_rff_3b, length_rff_4b, length_rff_5b, length_order,   &
                              ftnbody_omega_2b, ftnbody_omega_3b, ftnbody_omega_4b, ftnbody_omega_5b, &
                              ftnbody_phase_random_2b, ftnbody_phase_random_3b, ftnbody_phase_random_4b, & 
                              ftnbody_phase_random_5b, &
                              covar_matrix_2b, number_2b_ftnb, mean_2b, &
                              covar_matrix_3b, number_3b_ftnb, mean_3b, &
                              covar_matrix_4b, number_4b_ftnb, mean_4b, &
                              covar_matrix_5b, number_5b_ftnb, mean_5b, &
                              ltmp_mean_ftnbody, ftnbody_read_from_potential, &
                              ftnbody_model_id, FTNBODY_MODEL_GRAMM
    use module_kernel, only:  krff_type
    use module_sample_rand_ker, only: get_rff_diag_sigma
    use mld_logger
    use mld_mpi
    use mld_string, only: vtoa
    use my_mpi_subroutines, only: subworlds_allreduce_from_evrywhere_double, &
                                  subworlds_allreduce_from_evrywhere_vect_double, & 
                                  subworlds_allreduce_from_evrywhere_matrix_double
    implicit none
    integer :: iiii, ii 
    real(kind_double) :: mu
    real(kind_double), allocatable, dimension(:) :: total_mean 
    real(kind_double), allocatable, dimension(:,:) :: total_cov 
    character(len=80) :: chstring, string  

    _NAMECURRENT_("init_length_ftnbody_")
    _MLD_BEGIN_

    ! Prediction / MD from a saved potential: omega/phase come from random.xml
    ! (loaded in init_ftnbody); do not re-estimate the length scales here.
    if (ftnbody_read_from_potential) then
      call log_info("ML: ... ftnbody length estimation skipped (loaded from random.xml)")
      _MLD_END_
      return
    end if

    mu = 0.d0

    ! First step computing the means ...........................
    init_mode_ftnbody = .true.
    ltmp_mean_ftnbody = .true. 
    call log_info("ML: ... ftnbody begins , in "//NAMECURRENT//" estimation of length, step 1. the mean " )
    call main_compute_descriptors()
    if (l_body_order(2)) then
      number_2b_ftnb = 0.d0 
 
      call subworlds_allreduce_from_evrywhere_double(count_2b, number_2b_ftnb)
      
      if (allocated(total_mean)) deallocate(total_mean) 
      allocate(total_mean(size(mean_2b,1)))
      
      call subworlds_allreduce_from_evrywhere_vect_double(mean_2b, total_mean)

      mean_2b = total_mean / number_2b_ftnb 

      write(chstring, '(i10)') int(number_2b_ftnb)
      call log_info("ML: ... 2b number ftnbody  : "//trim(chstring))
      string=''
      do ii = 1, size(mean_2b,1)
        write(chstring, '(es20.9)') mean_2b(ii)
        string = trim(string)//trim(chstring)//" "
      end do   
      call log_info("ML: ... 2b mean   ftnbody  : "//trim(string))
      count_2b = 0 
      covar_matrix_2b(:,:)=0.d0 
    end if

    if (l_body_order(3)) then
      number_3b_ftnb = 0.d0 
      

      call subworlds_allreduce_from_evrywhere_double(count_3b, number_3b_ftnb)

      if (allocated(total_mean)) deallocate(total_mean) 
      allocate(total_mean(size(mean_3b,1)))


      call subworlds_allreduce_from_evrywhere_vect_double(mean_3b, total_mean)

      mean_3b = total_mean / number_3b_ftnb 

      write(chstring, '(i10)') int(number_3b_ftnb)
      call log_info("ML: ... 3b number ftnbody  : "//trim(chstring))

      string=''
      do ii = 1, size(mean_3b,1)
        write(chstring, '(es20.9)') mean_3b(ii)
        string = trim(string)//trim(chstring)//" "
      end do  
      call log_info("ML: ... 3b mean   ftnbody  : "//trim(string))

      count_3b = 0 
      covar_matrix_3b(:,:)=0.d0 
    end if 

    if (l_body_order(4)) then
      number_4b_ftnb = 0.d0 
      

      call subworlds_allreduce_from_evrywhere_double(count_4b, number_4b_ftnb)

      if (allocated(total_mean)) deallocate(total_mean) 
      allocate(total_mean(size(mean_4b,1)))

      call subworlds_allreduce_from_evrywhere_vect_double(mean_4b, total_mean)

      mean_4b = total_mean / number_4b_ftnb 

      call log_info("ML: ... 4b number ftnbody  : "//vtoa(number_4b_ftnb))
      call log_info("ML: ... 4b mean   ftnbody  : "//vtoa(mean_4b))


      count_4b = 0 
      covar_matrix_4b(:,:)=0.d0 
    end if 
    if (l_body_order(5)) then
      number_5b_ftnb = 0.d0 
      
      call subworlds_allreduce_from_evrywhere_double(count_5b, number_5b_ftnb)

      if (allocated(total_mean)) deallocate(total_mean) 
      allocate(total_mean(size(mean_5b,1)))


      call subworlds_allreduce_from_evrywhere_vect_double(mean_5b, total_mean)

      mean_5b = total_mean / number_5b_ftnb 

      call log_info("ML: ... 5b number ftnbody  : "//vtoa(number_5b_ftnb))
      call log_info("ML: ... 5b mean   ftnbody  : "//vtoa(mean_5b))

      count_5b = 0
      covar_matrix_5b(:,:)=0.d0 
    end if 

    !  step computing the covariance matrix  ...........................
    call log_info("ML: ... ftnbody mean estimated , in "//NAMECURRENT//" estimation of length, step 2. the covariance " )
    init_mode_ftnbody = .true.
    ltmp_mean_ftnbody = .false. 
    call main_compute_descriptors()    
    init_mode_ftnbody = .false. 
    ltmp_mean_ftnbody = .false. 

    call log_info("ML: ... ftnbody completed the covariance " )

    if (l_body_order(2)) then
      if (allocated(total_cov)) deallocate(total_cov) 
      allocate(total_cov(size(covar_matrix_2b,1), size(covar_matrix_2b,2)))

      call subworlds_allreduce_from_evrywhere_matrix_double(covar_matrix_2b, total_cov)

      covar_matrix_2b(:,:) = total_cov(:,:) / number_2b_ftnb
      do iiii = 1, dim_qbody(2)
        length_rff_2b(iiii) = SQRT(covar_matrix_2b(iiii,iiii))
      end do
      length_rff_2b(:) = length_rff_2b(:)/2



      string=''
      do ii = 1, size(length_rff_2b,1)
        write(chstring, '(es20.9)') length_rff_2b(ii)
        string = trim(string)//trim(chstring)//" "
      end do  
      call log_info("ML: ... 2b length   ftnbody  : "//trim(string))

      length_order = 2
      call  get_rff_diag_sigma(krff_type, mu, length_rff_2b, dim_qbody(2), dim_rff(2), ftnbody_omega_2b, ftnbody_phase_random_2b)
    end if

    if (l_body_order(3)) then
      if (allocated(total_cov)) deallocate(total_cov) 
      allocate(total_cov(size(covar_matrix_3b,1), size(covar_matrix_3b,2)))
 

      call subworlds_allreduce_from_evrywhere_matrix_double(covar_matrix_3b, total_cov)

      covar_matrix_3b(:,:) = total_cov(:,:) / number_3b_ftnb

      do iiii = 1, dim_qbody(3)
        length_rff_3b(iiii) = SQRT(covar_matrix_3b(iiii,iiii))
      end do
      length_rff_3b(:) = length_rff_3b(:)/2
      ! Gram model: symmetry-related slots share one Fourier length scale
      if (ftnbody_model_id == FTNBODY_MODEL_GRAMM) call pool_gram_lengths(2, length_rff_3b)

      string=''
      do ii = 1, size(length_rff_3b,1)
        write(chstring, '(es20.9)') length_rff_3b(ii)
        string = trim(string)//trim(chstring)//" "
      end do  
      call log_info("ML: ... 3b length   ftnbody  : "//trim(string))

      length_order = 3
      call  get_rff_diag_sigma(krff_type, mu, length_rff_3b, dim_qbody(3), dim_rff(3), ftnbody_omega_3b, ftnbody_phase_random_3b)
    end if 

    if (l_body_order(4)) then
      if (allocated(total_cov)) deallocate(total_cov) 
      allocate(total_cov(size(covar_matrix_4b,1), size(covar_matrix_4b,2)))

      call subworlds_allreduce_from_evrywhere_matrix_double(covar_matrix_4b, total_cov)

      covar_matrix_4b(:,:) = total_cov(:,:) / number_4b_ftnb

      !covar_matrix_4b(:,:) = covar_matrix_4b(:,:) / SQRT(count_4b)

      do iiii = 1, dim_qbody(4)
        length_rff_4b(iiii) = SQRT(covar_matrix_4b(iiii,iiii))
      end do
      length_rff_4b(:) = length_rff_4b(:)/2
      ! Gram model: symmetry-related slots share one Fourier length scale
      if (ftnbody_model_id == FTNBODY_MODEL_GRAMM) call pool_gram_lengths(3, length_rff_4b)

      call log_info("ML: ... 4b length   ftnbody  : "//vtoa(length_rff_4b))

      length_order = 4
      call  get_rff_diag_sigma(krff_type, mu, length_rff_4b, dim_qbody(4), dim_rff(4), ftnbody_omega_4b, ftnbody_phase_random_4b)
    end if

    if (l_body_order(5)) then
      if (allocated(total_cov)) deallocate(total_cov) 
      allocate(total_cov(size(covar_matrix_5b,1), size(covar_matrix_5b,2)))


      call subworlds_allreduce_from_evrywhere_matrix_double(covar_matrix_5b, total_cov)
      
      covar_matrix_5b(:,:) = total_cov(:,:) / number_5b_ftnb

      !covar_matrix_5b(:,:) = covar_matrix_5b(:,:) / SQRT(count_5b)

      do iiii = 1, dim_qbody(5)
        length_rff_5b(iiii) = SQRT(covar_matrix_5b(iiii,iiii))
      end do
      length_rff_5b(:) = length_rff_5b(:)/2
      ! Gram model: symmetry-related slots share one Fourier length scale
      if (ftnbody_model_id == FTNBODY_MODEL_GRAMM) call pool_gram_lengths(4, length_rff_5b)

      call log_info("ML: ... 5b length   ftnbody  : "//vtoa(length_rff_5b))


      length_order = 5
      call  get_rff_diag_sigma(krff_type, mu, length_rff_5b, dim_qbody(5), dim_rff(5), ftnbody_omega_5b, ftnbody_phase_random_5b)
    end if

    _MLD_END_

  contains

    ! Gram model (ftnbody_model='gramm'): the m radial slots u_1..u_m are
    ! symmetry-related, and so are the C(m,2) angular slots c_pq; use a single
    ! pooled Fourier length l_u for the radial group and l_c for the angular one
    ! (docs/MLT5 sec. "Random-frequency scales").
    subroutine pool_gram_lengths(m, lvec)
      integer, intent(in) :: m
      real(kind_double), dimension(:), intent(inout) :: lvec
      real(kind_double) :: lu, lc
      integer :: d
      d = size(lvec, 1)
      lu = sum(lvec(1:m))/dble(m)
      lvec(1:m) = lu
      if (d > m) then
        lc = sum(lvec(m+1:d))/dble(d - m)
        lvec(m+1:d) = lc
      end if
    end subroutine pool_gram_lengths

  end subroutine init_length_ftnbody



  module module_compute_body_order
    use module_kind_variables, ONLY: kind_double
    implicit none 
    contains 
    subroutine ftnbody_order_2(type_db_ja, max_neigh_local, i_type_db, i_central, r_central, tmp_dxp, desc_forces_local)
      use module_ftnbody, only: tmp_real, tmp_dreal, ftnbody_omega_2b, ftnbody_phase_random_2b, delta_rff, dim_rff, &
                                mean_2b, init_mode_ftnbody, covar_matrix_2b , count_2b, ltmp_mean_ftnbody, &
                                r_cut_ft2b, r_cut_width_ft2b, &
                                ftnbody_chem_mode, ftnbody_n_channels, ftnbody_chem_center, ftnbody_chem_neigh
      use module_chemical_species, only: fix_no_of_elements
      use module_body_desc, only : dim_desc_body
      use module_neigh_local, only: r_cut_in, r_cut_width_in, r_cut_pair_in, fcut_rij_inout
      !use module_kernel, only:  length_kernel, sigma_kernel
      implicit none
      integer, intent(in) ::  type_db_ja,  max_neigh_local
      integer, dimension(:), intent(in)  :: i_type_db, i_central
      real(kind_double), dimension(:), intent(in) :: r_central
      real(kind_double), dimension(:,:), intent(in) :: tmp_dxp
      logical, intent(in) :: desc_forces_local
      integer :: ia, ix, ii
      real(kind_double), dimension(:), allocatable :: tmpf
      real(kind_double) :: norm
      real(kind_double), dimension(:), allocatable :: fcut_all, dfcut_all
      real(kind_double) :: fcut, dfcut
      real(kind_double) :: local_r_cut_in
      integer :: type_fcut_in, type_fcut_out
      ! --- multispecies chemistry (docs/perspective_ftnbody.md) ---
      integer :: nrff, nsp, nchan, sa, sj, ic, off, irk
      real(kind_double) :: sgn
      integer(8) :: hkey
      real(kind_double), dimension(:), allocatable :: geo, dgeo, chem_w
      ! hoisted cos/sin of the Fourier phase (evaluated once per neighbour)
      real(kind_double), dimension(:), allocatable :: cosv, sinv
      !TODOftnbody...
      integer :: izozo 
      integer, dimension(size(i_type_db)) :: iv1_zozo  
      integer, dimension(size(i_central)) :: iv2_zozo  
 
      izozo =  type_db_ja
      iv1_zozo = i_type_db
      iv2_zozo = i_central


      nrff  = dim_rff(2)                 ! geometry random-Fourier features (per channel)
      nsp   = fix_no_of_elements
      nchan = ftnbody_n_channels(2)      ! 1 (mode 0), S*S (mode 1), R (mode 2), H (mode 3)

      allocate(tmpf(nrff))
      allocate(geo(nrff))
      allocate(dgeo(nrff))
      allocate(chem_w(max(nchan,1)))
      allocate(cosv(nrff))
      allocate(sinv(nrff))
      allocate(fcut_all(max_neigh_local))
      allocate(dfcut_all(max_neigh_local))

      !TODOftnbody ...
      type_fcut_in = 3
      type_fcut_out = 2
      ! normalization uses the geometry feature count (dim_rff), not the expanded block
      norm = sqrt(1.d0)/sqrt(dble(nrff))*delta_rff(2)
      sa = type_db_ja
      ! compute_qja_a_partir_de_rcentral
      tmp_real(:) = 0.d0
      ! modes 1/2/3 write only the active channel block per neighbour -> zero the rest
      if (ftnbody_chem_mode /= 0 .and. desc_forces_local) tmp_dreal(:, :, :) = 0.d0
      do ii = 1, max_neigh_local
        ! Pair-specific inner cutoff: i_type_db(ii+1) due to assumed-shape remapping (index 1 = central)
        local_r_cut_in = r_cut_pair_in(type_db_ja, i_type_db(ii + 1))
        call fcut_rij_inout(type_fcut_in, type_fcut_out, r_central(ii), r_cut_ft2b, r_cut_width_ft2b, local_r_cut_in, r_cut_width_in, desc_forces_local, fcut, dfcut)
        fcut_all(ii) = fcut
        dfcut_all(ii) = dfcut
      end do

      do ia = 1, max_neigh_local
          if (init_mode_ftnbody) then
            count_2b = count_2b + 1
            if (ltmp_mean_ftnbody) then
              mean_2b(:) = mean_2b(:) + r_central(ia)
            else
              covar_matrix_2b(:,:) = covar_matrix_2b(:,:) + (r_central(ia) - mean_2b(1) ) **2
            end if
            cycle
          end if
          tmpf(1:nrff) = r_central(ia)*ftnbody_omega_2b(1,1:nrff) + ftnbody_phase_random_2b(1:nrff)
          ! evaluate the transcendentals once per neighbour (bit-identical hoist)
          cosv(1:nrff) = cos(tmpf(1:nrff))
          if (desc_forces_local) sinv(1:nrff) = sin(tmpf(1:nrff))

          fcut = fcut_all(ia)
          dfcut = dfcut_all(ia)

          if (ftnbody_chem_mode == 0) then
            ! ---- single species: bit-identical to the original implementation ----
            tmp_real(1:nrff) = tmp_real(1:nrff) + cosv(1:nrff) * norm * fcut
            if (desc_forces_local) then
              do ix = 1, 3
                tmp_dreal(1:nrff, ia, ix) =  - sinv(1:nrff) * ftnbody_omega_2b(1,1:nrff) * tmp_dxp(ix, ia+1)/r_central(ia)*norm * fcut &
                                             + cosv(1:nrff) * tmp_dxp(ix, ia+1)/r_central(ia) * norm * dfcut
              end do
            end if
            cycle
          end if

          ! ---- multispecies chemistry modes (1,2,3): geometry feature + channel routing ----
          geo(1:nrff) = cosv(1:nrff) * norm * fcut
          sj = i_type_db(ia + 1)
          call ftnbody_chem_route_2b(ftnbody_chem_mode, sa, sj, nsp, nrff, nchan, ic, off, sgn, chem_w)

          select case (ftnbody_chem_mode)
          case (2)   ! low-rank: spread over ranks with chemical weights
            do irk = 1, nchan
              off = (irk - 1)*nrff
              tmp_real(off+1:off+nrff) = tmp_real(off+1:off+nrff) + chem_w(irk)*geo(1:nrff)
            end do
          case default   ! mode 1 (exact), mode 3 (hash): one channel, weight sgn
            tmp_real(off+1:off+nrff) = tmp_real(off+1:off+nrff) + sgn*geo(1:nrff)
          end select

          if (desc_forces_local) then
            do ix = 1, 3
              dgeo(1:nrff) = - sinv(1:nrff) * ftnbody_omega_2b(1,1:nrff) * tmp_dxp(ix, ia+1)/r_central(ia)*norm * fcut &
                             + cosv(1:nrff) * tmp_dxp(ix, ia+1)/r_central(ia) * norm * dfcut
              select case (ftnbody_chem_mode)
              case (2)
                do irk = 1, nchan
                  off = (irk - 1)*nrff
                  tmp_dreal(off+1:off+nrff, ia, ix) = chem_w(irk)*dgeo(1:nrff)
                end do
              case default
                tmp_dreal(off+1:off+nrff, ia, ix) = sgn*dgeo(1:nrff)
              end select
            end do
          end if
      end do

      deallocate(tmpf, geo, dgeo, chem_w, cosv, sinv)
    end subroutine ftnbody_order_2

    ! ------------------------------------------------------------------
    ! Resolve the chemical channel routing of one 2-body cluster (central
    ! species sa, neighbour species sj) for each ftnbody_chem_mode.
    ! Returns: ic (channel id, 1-based), off=(ic-1)*nrff, sgn (channel weight),
    ! and chem_w(1:nchan) for the low-rank mode.  See docs/perspective_ftnbody.md.
    ! ------------------------------------------------------------------
    subroutine ftnbody_chem_route_2b(mode, sa, sj, nsp, nrff, nchan, ic, off, sgn, chem_w)
      use module_ftnbody, only: ftnbody_chem_center, ftnbody_chem_neigh
      implicit none
      integer, intent(in)  :: mode, sa, sj, nsp, nrff, nchan
      integer, intent(out) :: ic, off
      real(kind_double), intent(out) :: sgn
      real(kind_double), dimension(:), intent(inout) :: chem_w
      integer :: irk
      integer(8) :: hkey
      sgn = 1.d0
      ic  = 1
      select case (mode)
      case (0)               ! single species
        ic = 1
      case (1)               ! exact channel (sa, sj), 1..nsp*nsp
        ic = (sa - 1)*nsp + sj
      case (2)               ! low-rank embedding weights (channel = rank block)
        do irk = 1, nchan
          chem_w(irk) = ftnbody_chem_center(sa, irk, 2) * ftnbody_chem_neigh(sj, irk, 1, 2)
        end do
        ic = 1
      case (3)               ! feature hashing of the (sa,sj) tuple into nchan channels
        hkey = int(sa - 1, 8)*int(nsp, 8) + int(sj - 1, 8)
        ic  = int(mod(hkey*2654435761_8, int(nchan, 8)), kind=4) + 1
        if (mod(hkey*40503_8 + 12345_8, 2_8) == 0_8) then
          sgn = 1.d0
        else
          sgn = -1.d0
        end if
      end select
      off = (ic - 1)*nrff
    end subroutine ftnbody_chem_route_2b


subroutine ftnbody_order_3(type_db_ja, max_neigh_local, i_type_db, i_central, r_central, tmp_dxp, ur_central, d_ur_central, desc_forces_local)
    ! compute all the body triangle with the following structure. Please pay attention that the distances
     ! are already r_cut -ed.
     !
     !     i2
     !    /
     !  x1
     ! / x3
     ! j--x2--i3
      use module_ftnbody, only: tmp_real, tmp_dreal, ftnbody_omega_3b, ftnbody_phase_random_3b, delta_rff, dim_rff, &
                                init_mode_ftnbody, covar_matrix_3b , count_3b, dim_qbody, mean_3b, ltmp_mean_ftnbody, &
                                r_cut_ft3b, r_cut_width_ft3b, &
                                ftnbody_chem_mode, ftnbody_n_channels, ftnbody_chem_center, ftnbody_chem_neigh
      use module_chemical_species, only: fix_no_of_elements
      use module_body_desc, only : dim_desc_body
      use module_neigh_local, only: r_cut_in, r_cut_width_in, r_cut_pair_in, fcut_rij_inout
      !use module_kernel, only:  length_kernel, sigma_kernel
      implicit none
      integer, intent(in) ::  type_db_ja,  max_neigh_local
      integer, dimension(:), intent(in)  :: i_type_db, i_central
      real(kind_double), dimension(:), intent(in) :: r_central, ur_central
      real(kind_double), dimension(:,:), intent(in) :: tmp_dxp
      real(kind_double), dimension(:,:), intent(in) :: d_ur_central
      logical, intent(in) :: desc_forces_local
      integer :: i2, i3, ix, kk, ii, dimf
      integer :: iiii, jjjj
      real(kind_double), dimension(:), allocatable :: tmpf
      integer, parameter :: dimq = 3
      real(kind_double), dimension(dimq)   :: qsym
      real(kind_double), dimension(3, dimq)      :: d2_qsym, d3_qsym
      real(kind_double), dimension(dim_rff(3),dimq) :: sin_tmp
      real(kind_double), dimension(dim_rff(3)) :: d2t, d3t
      ! hoisted cos/sin of the Fourier phase (evaluated once per cluster)
      real(kind_double), dimension(dim_rff(3)) :: cosv, sinv
      ! --- multispecies chemistry (docs/perspective_ftnbody.md) ---
      integer :: nrff, nsp, nchan, sa, sj1, sj2, ic, off, irk
      real(kind_double) :: sgn
      real(kind_double), dimension(dim_rff(3)) :: gE, gF2, gF3
      real(kind_double), dimension(:), allocatable :: chem_w
      real(kind_double) :: cos_2j3, x1, x2, x3, norm, one=1.d0 !, one_to_r_central_i2=0.d0
      integer :: incx = 1, incy = 1 
      real(kind_double), dimension(3)   :: d2_x1, d2_x2, d2_x3, &
                                           d3_x1, d3_x2, d3_x3
      real(kind_double), dimension(:), allocatable :: fcut_all, dfcut_all 
      real(kind_double) :: fcut, dfcut, fcut_i2, dfcut_i2, fcut_i3, dfcut_i3
      real(kind_double) :: local_r_cut_in
      integer :: type_fcut_in, type_fcut_out 
      !TODOftnbody...
      integer :: izozo 
      integer, dimension(size(i_type_db)) :: iv1_zozo  
      integer, dimension(size(i_central)) :: iv2_zozo  
 
      izozo =  type_db_ja
      iv1_zozo = i_type_db 
      iv2_zozo = i_central


      
      nrff  = dim_rff(3)               ! geometry random-Fourier features (per channel)
      nsp   = fix_no_of_elements
      nchan = ftnbody_n_channels(3)    ! 1 (mode 0), N_chem (mode 1), R (mode 2), H (mode 3)
      sa    = type_db_ja

      allocate(tmpf(nrff))
      allocate(chem_w(max(nchan,1)))
      allocate(fcut_all(max_neigh_local))
      allocate(dfcut_all(max_neigh_local))
      ! allocate(ttmp(dim_desc_body(3),3))

      !TODOftnbody
      type_fcut_in=3
      type_fcut_out=2

      ! normalization uses the geometry feature count (dim_rff), not the expanded block
      norm = sqrt(1.d0)/sqrt(dble(nrff))*delta_rff(3)
      ! compute_qja_a_partir_de_rcentral
      dimf = nrff
      tmpf(:) = ftnbody_phase_random_3b

      if (max_neigh_local <= 1) return
      tmp_real(:) = 0.d0 
      if (desc_forces_local) then
        tmp_dreal(:, :, :) = 0.d0
    
        d2_x1(:) = 0.d0
        d2_x2(:) = 0.d0
        d2_x3(:) = 0.d0
    
        d3_x1(:) = 0.d0
        d3_x2(:) = 0.d0
        d3_x3(:) = 0.d0

      end if

      do ii = 1, max_neigh_local
        ! Pair-specific inner cutoff: i_type_db(ii+1) due to assumed-shape remapping (index 1 = central)
        local_r_cut_in = r_cut_pair_in(type_db_ja, i_type_db(ii + 1))
        call fcut_rij_inout(type_fcut_in, type_fcut_out, r_central(ii), r_cut_ft3b, r_cut_width_ft3b, local_r_cut_in, r_cut_width_in, desc_forces_local, fcut, dfcut)
        fcut_all(ii) = fcut
        dfcut_all(ii) = dfcut
      end do

      do i2 = 1, max_neigh_local

#ifdef MLD_NDM
        ! ttmp(:,:)= 0.d0 
        !one_to_r_central_i2 = 1.d0/r_central(i2)
#endif
        fcut_i2 = fcut_all(i2)
        if (fcut_i2 == 0.d0) cycle
        dfcut_i2 = dfcut_all(i2)

        x1 = ur_central(i2)
        if (desc_forces_local) then
          d2_x1(1:3) = d_ur_central(1:3, i2)
        end if
        do i3 = i2 + 1, max_neigh_local
            fcut_i3 = fcut_all(i3)
            if (fcut_i3 == 0.d0) cycle
            dfcut_i3 = dfcut_all(i3)
          ! if (i2/=i3) then  
            x2 = ur_central(i3)
            cos_2j3 = dot_product(tmp_dxp(1:3, i2+1), tmp_dxp(1:3, i3+1))/(r_central(i2)*r_central(i3))
            x3 = cos_2j3
            if (desc_forces_local) then
              d3_x2(1:3) = d_ur_central(1:3, i3)
              d2_x3(1:3) = tmp_dxp(1:3, i3+1)/(r_central(i2)*r_central(i3)) - tmp_dxp(1:3, i2+1)*cos_2j3/r_central(i2)**2
              d3_x3(1:3) = tmp_dxp(1:3, i2+1)/(r_central(i2)*r_central(i3)) - tmp_dxp(1:3, i3+1)*cos_2j3/r_central(i3)**2
            end if

            qsym(1) = x1 + x2
            qsym(2) = x1 * x2
            qsym(3) = x3 
           
            if (init_mode_ftnbody) then
              count_3b = count_3b + 1
              if (ltmp_mean_ftnbody) then 
                mean_3b(:) = mean_3b(:)  + qsym(:) 
              else 
                do iiii = 1, dim_qbody(3)
                  do jjjj = 1, dim_qbody(3)
                    covar_matrix_3b(iiii,jjjj) = covar_matrix_3b(iiii,jjjj) + ( qsym(iiii)  & 
                                                 - mean_3b(iiii))* (qsym(jjjj) - mean_3b(jjjj)) 
                  end do
                end do
              end if 
              cycle
            end if

            if (desc_forces_local) then
              ! 2 is in 1 - dist ;  3 - angles
              d2_qsym(1:3, 1) = d2_x1(1:3)
              d2_qsym(1:3, 2) = x2 * d2_x1(1:3)
              d2_qsym(1:3, 3) = d2_x3(1:3)
    
              ! 3 is in 2 - dist ;  3 - angles
              d3_qsym(1:3, 1) = d3_x2(1:3)
              d3_qsym(1:3, 2) = x1 * d3_x2(1:3)
              d3_qsym(1:3, 3) = d3_x3(1:3)
            end if

            tmpf(:) = ftnbody_phase_random_3b(:)
            call dgemv('T', dimq, dimf, one, ftnbody_omega_3b, dimq, qsym, incx, one, tmpf, incy)
            !!! tmpf(:) = qsym(1:3)*ftnbody_omega_3b(1:3,:) + ftnbody_phase_random_3b(:)

            ! evaluate the transcendentals once per cluster (bit-identical hoist);
            ! sin/sin_tmp are only needed for the force block
            cosv(:) = cos(tmpf(:))
            if (desc_forces_local) then
              sinv(:) = sin(tmpf(:))
              do ii = 1,dimq
                sin_tmp(:,ii) =   -sinv(:) * ftnbody_omega_3b(ii,:)
              end do
            end if

            if (ftnbody_chem_mode == 0) then
              ! ---- single species: bit-identical to the original implementation ----
              tmp_real(:) = tmp_real(:) + cosv(:) * fcut_i2 * fcut_i3
              if (desc_forces_local) then
                do ix =1,3
                  d2t(:) =                    sin_tmp(:,1) * d2_qsym(ix,1) &
                                            + sin_tmp(:,2) * d2_qsym(ix,2) &
                                            + sin_tmp(:,3) * d2_qsym(ix,3)

                  d3t(:) =                    sin_tmp(:,1) * d3_qsym(ix,1) &
                                            + sin_tmp(:,2) * d3_qsym(ix,2) &
                                            + sin_tmp(:,3) * d3_qsym(ix,3)

                  tmp_dreal(:,i2,ix) = tmp_dreal(:,i2,ix) + d2t(:) * fcut_i2 * fcut_i3 + cosv(:) * fcut_i3 * dfcut_i2 * tmp_dxp(ix, i2+1)/r_central(i2)
                  tmp_dreal(:,i3,ix) = tmp_dreal(:,i3,ix) + d3t(:) * fcut_i2 * fcut_i3 + cosv(:) * fcut_i2 * dfcut_i3 * tmp_dxp(ix, i3+1)/r_central(i3)

                end do
              end if
            else
              ! ---- multispecies modes 1/2/3: triangle (sa, sj1, sj2) channel routing ----
              sj1 = i_type_db(i2 + 1)
              sj2 = i_type_db(i3 + 1)
              call ftnbody_chem_route_3b(ftnbody_chem_mode, sa, sj1, sj2, nsp, nrff, nchan, ic, off, sgn, chem_w)
              gE(1:nrff) = cosv(1:nrff) * fcut_i2 * fcut_i3
              if (ftnbody_chem_mode == 2) then
                do irk = 1, nchan
                  off = (irk - 1)*nrff
                  tmp_real(off+1:off+nrff) = tmp_real(off+1:off+nrff) + chem_w(irk)*gE(1:nrff)
                end do
              else
                tmp_real(off+1:off+nrff) = tmp_real(off+1:off+nrff) + sgn*gE(1:nrff)
              end if
              if (desc_forces_local) then
                do ix = 1, 3
                  d2t(:) = sin_tmp(:,1)*d2_qsym(ix,1) + sin_tmp(:,2)*d2_qsym(ix,2) + sin_tmp(:,3)*d2_qsym(ix,3)
                  d3t(:) = sin_tmp(:,1)*d3_qsym(ix,1) + sin_tmp(:,2)*d3_qsym(ix,2) + sin_tmp(:,3)*d3_qsym(ix,3)
                  gF2(1:nrff) = d2t(1:nrff)*fcut_i2*fcut_i3 + cosv(1:nrff)*fcut_i3*dfcut_i2*tmp_dxp(ix, i2+1)/r_central(i2)
                  gF3(1:nrff) = d3t(1:nrff)*fcut_i2*fcut_i3 + cosv(1:nrff)*fcut_i2*dfcut_i3*tmp_dxp(ix, i3+1)/r_central(i3)
                  if (ftnbody_chem_mode == 2) then
                    do irk = 1, nchan
                      off = (irk - 1)*nrff
                      tmp_dreal(off+1:off+nrff,i2,ix) = tmp_dreal(off+1:off+nrff,i2,ix) + chem_w(irk)*gF2(1:nrff)
                      tmp_dreal(off+1:off+nrff,i3,ix) = tmp_dreal(off+1:off+nrff,i3,ix) + chem_w(irk)*gF3(1:nrff)
                    end do
                  else
                    tmp_dreal(off+1:off+nrff,i2,ix) = tmp_dreal(off+1:off+nrff,i2,ix) + sgn*gF2(1:nrff)
                    tmp_dreal(off+1:off+nrff,i3,ix) = tmp_dreal(off+1:off+nrff,i3,ix) + sgn*gF3(1:nrff)
                  end if
                end do
              end if
            end if
          ! end if
        end do
        ! do ix =1,3
        !    tmp_dreal(:,i2,ix) =  ttmp(:,ix)*norm
        ! end do 
      end do
      ! norm = norm!/2.d0      
      tmp_real(:) = tmp_real(:)*norm
      if (desc_forces_local) then 
        do kk = 1, max_neigh_local
          do ix =1, 3
            tmp_dreal(:,kk,ix) = tmp_dreal(:,kk,ix)*norm  
          end do       
        end do   ! end kk 
      end if 

      deallocate(tmpf, chem_w)
     end subroutine ftnbody_order_3

    ! ------------------------------------------------------------------
    ! Channel routing of one 3-body triangle (central species sa, neighbour
    ! species sj1, sj2) for each ftnbody_chem_mode.  The compact 3-body q is
    ! already symmetric in the two neighbours, so exact/hash channels use the
    ! sorted neighbour pair and the low-rank factor is symmetrized over the two
    ! neighbour permutations.  See docs/perspective_ftnbody.md.
    ! ------------------------------------------------------------------
    subroutine ftnbody_chem_route_3b(mode, sa, sj1, sj2, nsp, nrff, nchan, ic, off, sgn, chem_w)
      use module_ftnbody, only: ftnbody_chem_center, ftnbody_chem_neigh
      implicit none
      integer, intent(in)  :: mode, sa, sj1, sj2, nsp, nrff, nchan
      integer, intent(out) :: ic, off
      real(kind_double), intent(out) :: sgn
      real(kind_double), dimension(:), intent(inout) :: chem_w
      integer :: irk, p, q, mset, mpair
      integer(8) :: hkey
      sgn = 1.d0
      ic  = 1
      p = min(sj1, sj2)        ! unordered neighbour pair (p <= q)
      q = max(sj1, sj2)
      select case (mode)
      case (1)                 ! exact channel: central x unordered neighbour multiset
        mpair = nsp*(nsp + 1)/2                          ! # unordered pairs with repetition
        mset  = (p - 1)*(2*nsp - p + 2)/2 + (q - p + 1)  ! rank of {p,q} in 1..mpair
        ic = (sa - 1)*mpair + mset
      case (2)                 ! low-rank, symmetrized over the 2 neighbour permutations
        do irk = 1, nchan
          chem_w(irk) = ftnbody_chem_center(sa, irk, 3) * &
               ( ftnbody_chem_neigh(sj1, irk, 1, 3)*ftnbody_chem_neigh(sj2, irk, 2, 3) &
               + ftnbody_chem_neigh(sj2, irk, 1, 3)*ftnbody_chem_neigh(sj1, irk, 2, 3) )
        end do
        ic = 1
      case (3)                 ! hash of (sa, sorted pair)
        hkey = (int(sa - 1, 8)*int(nsp, 8) + int(p - 1, 8))*int(nsp, 8) + int(q - 1, 8)
        ic  = int(mod(hkey*2654435761_8, int(nchan, 8)), kind=4) + 1
        if (mod(hkey*40503_8 + 12345_8, 2_8) == 0_8) then
          sgn = 1.d0
        else
          sgn = -1.d0
        end if
      end select
      off = (ic - 1)*nrff
    end subroutine ftnbody_chem_route_3b

     subroutine ftnbody_order_4(type_db_ja, max_neigh_local, i_type_db, i_central, r_central, tmp_dxp, ur_central, d_ur_central, desc_forces_local)
     
     ! compute all the body4 triangle with the following structure. Please pay attention that the distances
     ! are already r_cut -ed.
     !
     !     i4
     !    /
     !  x3
     ! / x5
     ! j--x1--i2  x6
     ! \ x4
     !  x2
     !   \
     !    i3


      use module_ftnbody, only: tmp_real, tmp_dreal, ftnbody_omega_4b, ftnbody_phase_random_4b, delta_rff, dim_rff, &
                                init_mode_ftnbody, covar_matrix_4b , count_4b, dim_qbody, mean_4b, ltmp_mean_ftnbody, &
                                r_cut_ft4b, r_cut_width_ft4b, &
                                ftnbody_chem_mode, ftnbody_n_channels, ftnbody_chem_center, ftnbody_chem_neigh
      use module_chemical_species, only: fix_no_of_elements
      use time_check_general, only : debug_time, MY_MPI_WTIME
      use module_body_desc, only : dim_desc_body
      use module_ftnbody,   only : t4b_inner_init, t4b_inner_desc, t4b_inner_deriv, t4b_inner_last
      use module_neigh_local, only: r_cut_in, r_cut_width_in, r_cut_pair_in, fcut_rij_inout
      implicit none
      integer, intent(in) ::  type_db_ja,  max_neigh_local
      integer, dimension(:), intent(in)  :: i_type_db, i_central
      real(kind_double), dimension(:), intent(in) :: r_central, ur_central 
      real(kind_double), dimension(:,:), intent(in) :: tmp_dxp
      real(kind_double), dimension(:,:), intent(in) :: d_ur_central
      logical, intent(in) :: desc_forces_local 
      integer :: i2, i3, i4, ix, kk, ii,  dimf
      integer :: iiii, jjjj
      real(kind_double), dimension(:), allocatable :: tmpf
      ! real(kind_double), dimension(:,:), allocatable :: ttmp 
      real(kind_double) :: norm, one=1.d0 ,  time00, time11, time22, time33 
      integer :: incx = 1, incy = 1!, icount_nn
      !cos! real(kind_double), dimension(dim_desc_body(4), 3, max_neigh_local**3)    :: d2_tmpr, d3_tmpr, d4_tmpr
      ! integer, dimension(max_neigh_local, max_neigh_local, max_neigh_local) :: innvec
      real(kind_double), dimension(3)   :: d2_x1, d2_x2, d2_x3, d2_x4, d2_x5, d2_x6, &
                                           d3_x1, d3_x2, d3_x3, d3_x4, d3_x5, d3_x6, &
                                           d4_x1, d4_x2, d4_x3, d4_x4, d4_x5, d4_x6
      real(kind_double)        :: cos_2j3, cos_3j4, cos_2j4, x1, x2, x3, x4, x5, x6, &
                                  x1_2, x1_3, x2_2, x2_3, x3_2, x3_3, &
                                  x4_2, x4_3, x5_2, x5_3, x6_2, x6_3
      integer, parameter :: dimq =6
      real(kind_double), dimension(dimq)   :: qsym
      real(kind_double), dimension(3, dimq)      :: d2_qsym, d3_qsym, d4_qsym
      real(kind_double), dimension(dim_rff(4)) :: d2t, d3t, d4t
      ! hoisted cos/sin of the Fourier phase (evaluated once per cluster)
      real(kind_double), dimension(dim_rff(4)) :: cosv, sinv
      ! R3: slot contraction as one dgemm  Wmat = omega^T * DQmat,
      ! DQmat(i, (t-1)*3+ix) = d{t}_qsym(ix, i); -sin applied afterwards
      real(kind_double), dimension(dimq, 9) :: DQmat
      real(kind_double), dimension(dim_rff(4), 9) :: Wmat
      real(kind_double) :: zero = 0.d0
      ! --- multispecies chemistry (docs/perspective_ftnbody.md) ---
      integer :: nrff, nsp, nchan, sa, sj1, sj2, sj3, ic, off, irk
      real(kind_double) :: sgn
      real(kind_double), dimension(dim_rff(4)) :: gE, gF2, gF3, gF4
      real(kind_double), dimension(:), allocatable :: chem_w
      real(kind_double), dimension(:), allocatable :: fcut_all, dfcut_all
      real(kind_double) :: fcut, dfcut, fcut_i2, dfcut_i2, fcut_i3, dfcut_i3, fcut_i4, dfcut_i4
      real(kind_double) :: local_r_cut_in
      integer :: type_fcut_in, type_fcut_out 
      !TODOftnbody...
      integer :: izozo 
      integer, dimension(size(i_type_db)) :: iv1_zozo  
      integer, dimension(size(i_central)) :: iv2_zozo  
 
       izozo =  type_db_ja
       iv1_zozo = i_type_db
       iv2_zozo = i_central


      !TODOftnbody ... 
      type_fcut_in=3
      type_fcut_out=2 

      nrff  = dim_rff(4)
      nsp   = fix_no_of_elements
      nchan = ftnbody_n_channels(4)
      sa    = type_db_ja

      allocate(tmpf(nrff))
      allocate(chem_w(max(nchan,1)))
      ! allocate(ttmp(dim_desc_body(4),dimq))
      allocate(fcut_all(max_neigh_local))
      allocate(dfcut_all(max_neigh_local))

      ! normalization uses the geometry feature count (dim_rff), not the expanded block
      norm = sqrt(1.d0)/sqrt(dble(nrff))*delta_rff(4)
      ! compute_qja_a_partir_de_rcentral
      dimf = nrff
      tmpf(:) = ftnbody_phase_random_4b

      if (max_neigh_local <= 1) return
      tmp_real(:) = 0.d0
      if (desc_forces_local) then
        !cos! d2_tmpr(:, :, :) = 0.d0
        !cos! d3_tmpr(:, :, :) = 0.d0
        !cos! d4_tmpr(:, :, :) = 0.d0
        tmp_dreal(:, :, :) = 0.d0
    
        d2_x1(:) = 0.d0
        d2_x2(:) = 0.d0
        d2_x3(:) = 0.d0
        d2_x4(:) = 0.d0
        d2_x5(:) = 0.d0
        d2_x6(:) = 0.d0
    
    
        d3_x1(:) = 0.d0
        d3_x2(:) = 0.d0
        d3_x3(:) = 0.d0
        d3_x4(:) = 0.d0
        d3_x5(:) = 0.d0
        d3_x6(:) = 0.d0
    
        d4_x1(:) = 0.d0
        d4_x2(:) = 0.d0
        d4_x3(:) = 0.d0
        d4_x4(:) = 0.d0
        d4_x5(:) = 0.d0
        d4_x6(:) = 0.d0
      end if

      do ii = 1, max_neigh_local
        ! Pair-specific inner cutoff: i_type_db(ii+1) due to assumed-shape remapping (index 1 = central)
        local_r_cut_in = r_cut_pair_in(type_db_ja, i_type_db(ii + 1))
        call fcut_rij_inout(type_fcut_in, type_fcut_out, r_central(ii), r_cut_ft4b, r_cut_width_ft4b, local_r_cut_in, r_cut_width_in, desc_forces_local, fcut, dfcut)
        fcut_all(ii) = fcut
        dfcut_all(ii) = dfcut
      end do

      ! icount_nn = 0 
      do i2 = 1, max_neigh_local
        fcut_i2 = fcut_all(i2)
        if (fcut_i2 == 0.d0) cycle
        dfcut_i2 = dfcut_all(i2)

        x1 = ur_central(i2)
        x1_2 = x1**2
        x1_3 = x1_2*x1
        if (desc_forces_local) then
          d2_x1(1:3) = d_ur_central(1:3, i2)
        end if
        do i3 = i2 + 1, max_neigh_local
          fcut_i3 = fcut_all(i3)
          if (fcut_i3 == 0.d0) cycle
          dfcut_i3 = dfcut_all(i3)
          
          x2 = ur_central(i3)
          x2_2 = x2**2
          x2_3 = x2_2*x2
          cos_2j3 = dot_product(tmp_dxp(:, i2+1), tmp_dxp(:, i3+1))/(r_central(i2)*r_central(i3))
          x4 = cos_2j3
          x4_2 = x4**2
          x4_3 = x4_2*x4
          if (desc_forces_local) then
            d3_x2(1:3) = d_ur_central(1:3, i3)
            d2_x4(1:3) = tmp_dxp(1:3, i3+1)/(r_central(i2)*r_central(i3)) - tmp_dxp(1:3, i2+1)*cos_2j3/r_central(i2)**2
            d3_x4(1:3) = tmp_dxp(1:3, i2+1)/(r_central(i2)*r_central(i3)) - tmp_dxp(1:3, i3+1)*cos_2j3/r_central(i3)**2
          end if
    
          ! ttmp(:,:)= 0.d0 
          do i4 = i3 + 1, max_neigh_local
            fcut_i4 = fcut_all(i4)
            if (fcut_i4 == 0.d0) cycle
            dfcut_i4 = dfcut_all(i4)

            if (debug_time) time00 = MY_MPI_WTIME()
            x3 = ur_central(i4)
            x3_2 = x3**2
            x3_3 = x3_2*x3
            cos_3j4 = dot_product(tmp_dxp(:, i3+1), tmp_dxp(:, i4+1))/(r_central(i3)*r_central(i4))
            x6 = cos_3j4
            x6_2 = x6**2
            x6_3 = x6_2*x6    
            cos_2j4 = dot_product(tmp_dxp(:, i2+1), tmp_dxp(:, i4+1))/(r_central(i2)*r_central(i4))
            x5 = cos_2j4
            x5_2 = x5**2
            x5_3 = x5_2*x5

            if (desc_forces_local) then
              d4_x3(1:3) = d_ur_central(1:3, i4)
              d3_x6(1:3) = tmp_dxp(1:3, i4+1)/(r_central(i3)*r_central(i4)) - tmp_dxp(1:3, i3+1)*cos_3j4/r_central(i3)**2
              d4_x6(1:3) = tmp_dxp(1:3, i3+1)/(r_central(i3)*r_central(i4)) - tmp_dxp(1:3, i4+1)*cos_3j4/r_central(i4)**2
              d2_x5(1:3) = tmp_dxp(1:3, i4+1)/(r_central(i2)*r_central(i4)) - tmp_dxp(1:3, i2+1)*cos_2j4/r_central(i2)**2
              d4_x5(1:3) = tmp_dxp(1:3, i2+1)/(r_central(i2)*r_central(i4)) - tmp_dxp(1:3, i4+1)*cos_2j4/r_central(i4)**2
            end if
            ! AAA case
            ! if ( (i_type(i2)==i_type(i3)).and.(i_type(i2)==i_type(i4)) ) then
            ! qsym(1) = x1 + x2 + x3
            ! qsym(2) = x4 + x5 + x6
            ! !qsym(3) = x1**2 + x2**2 + x3**2
            ! qsym(3) = x1_2 + x2_2 + x3_2
            ! ! qsym(4) = x1*x4 + x2*x5 + x3*x6
            ! qsym(4) = qsym(1) * qsym(2)!(x1 + x2 + x3) * (x4 + x5 + x6)
            ! qsym(5) = x4**3 + x5**3 + x6**3
            ! ! qsym(6) = x1**3 + x2**3 + x3**3 + x4**2*x5 + x4*x6**2 + x6*x5**2
            ! ! qsym(6) = x1_3 + x2_3 + x3_3 + x4**2*x5 + x4*x6**2 + x6*x5**2
            ! qsym(6) = x1_3 + x2_3 + x3_3 + x4*x5*x6

            qsym(1) = x1 + x2 + x3
            qsym(2) = x4 + x5 + x6
            qsym(3) = x1_2 + x2_2 + x3_2
            qsym(4) = x4_2 + x5_2 + x6_2
            qsym(5) = x1_3 + x2_3 + x3_3  
            qsym(6) = x4_3 + x5_3 + x6_3
            
            if (init_mode_ftnbody) then
              count_4b = count_4b + 1
              if (ltmp_mean_ftnbody) then 
                mean_4b(:) = mean_4b(:) + qsym(:)  
              else 
                do iiii = 1, dim_qbody(4)
                  do jjjj = 1, dim_qbody(4) 
                    covar_matrix_4b(iiii,jjjj) = covar_matrix_4b(iiii,jjjj) + (qsym(iiii) - mean_4b(iiii)) &
                                                                             *(qsym(jjjj) - mean_4b(jjjj))
                  end do
                end do
              end if 
              cycle
            end if

            if (desc_forces_local) then
    
              ! 2 is in 1 - dist ;  4,5 - angles
              d2_qsym(1:3, 1) = d2_x1(1:3)
              d2_qsym(1:3, 2) = d2_x4(1:3) + d2_x5(1:3)
              d2_qsym(1:3, 3) = 2.d0*x1*d2_x1(1:3)
              d2_qsym(1:3, 4) = 2.d0*x4*d2_x4(1:3) + 2.d0*x5*d2_x5(1:3)
              d2_qsym(1:3, 5) = 3.d0*x1_2*d2_x1(1:3)
              d2_qsym(1:3, 6) = 3.d0*x4_2*d2_x4(1:3) + 3.d0*x5_2*d2_x5(1:3)
    
              ! 3 is in 2 - dist ;  4,6 - angles
              d3_qsym(1:3, 1) = d3_x2(1:3)
              d3_qsym(1:3, 2) = d3_x4(1:3) + d3_x6(1:3)
              d3_qsym(1:3, 3) = 2.d0*x2*d3_x2(1:3)
              d3_qsym(1:3, 4) = 2.d0*x4*d3_x4(1:3) + 2.d0*x6*d3_x6(1:3)
              d3_qsym(1:3, 5) = 3.d0*x2_2*d3_x2(1:3)
              d3_qsym(1:3, 6) = 3.d0*x4_2*d3_x4(1:3) + 3.d0*x6_2*d3_x6(1:3)

              ! 4 is in 3 - dist ;  5,6 - angles
              d4_qsym(1:3, 1) = d4_x3(1:3)
              d4_qsym(1:3, 2) = d4_x5(1:3) + d4_x6(1:3)
              d4_qsym(1:3, 3) = 2.d0*x3*d4_x3(1:3)
              d4_qsym(1:3, 4) = 2.d0*x5*d4_x5(1:3) + 2.d0*x6*d4_x6(1:3)
              d4_qsym(1:3, 5) = 3.d0*x3_2*d4_x3(1:3)
              d4_qsym(1:3, 6) = 3.d0*x5_2*d4_x5(1:3) + 3.d0*x6_2*d4_x6(1:3)
              
            end if
            if (debug_time) time11 = MY_MPI_WTIME()

            tmpf(:) = ftnbody_phase_random_4b(:)
            call dgemv('T', dimq, dimf, one, ftnbody_omega_4b, dimq, qsym, incx, one, tmpf, incy)

            if (debug_time) time22 = MY_MPI_WTIME()

            ! evaluate the transcendentals once per cluster; R3: the slot
            ! contraction is one dgemm on the raw omega, -sin applied after
            cosv(:) = cos(tmpf(:))
            if (desc_forces_local) then
              sinv(:) = sin(tmpf(:))
              do ii = 1, dimq
                DQmat(ii, 1:3) = d2_qsym(1:3, ii)
                DQmat(ii, 4:6) = d3_qsym(1:3, ii)
                DQmat(ii, 7:9) = d4_qsym(1:3, ii)
              end do
              call dgemm('T', 'N', nrff, 9, dimq, one, ftnbody_omega_4b, dimq, DQmat, dimq, zero, Wmat, nrff)
            end if

            if (ftnbody_chem_mode == 0) then
              ! ---- single species: bit-identical to the original implementation ----
              tmp_real(:) = tmp_real(:) + cosv(:) * fcut_i2 * fcut_i3 * fcut_i4
              if (desc_forces_local) then
                do ix = 1, 3
                  d2t(:) = -sinv(:)*Wmat(:,ix)
                  d3t(:) = -sinv(:)*Wmat(:,3+ix)
                  d4t(:) = -sinv(:)*Wmat(:,6+ix)

                  tmp_dreal(:,i2,ix) = tmp_dreal(:,i2,ix) + d2t(:) * fcut_i2 * fcut_i3 * fcut_i4 + cosv(:) * dfcut_i2 * fcut_i3 * fcut_i4 * tmp_dxp(ix, i2+1)/r_central(i2)
                  tmp_dreal(:,i3,ix) = tmp_dreal(:,i3,ix) + d3t(:) * fcut_i2 * fcut_i3 * fcut_i4 + cosv(:) * fcut_i2 * dfcut_i3 * fcut_i4 * tmp_dxp(ix, i3+1)/r_central(i3)
                  tmp_dreal(:,i4,ix) = tmp_dreal(:,i4,ix) + d4t(:) * fcut_i2 * fcut_i3 * fcut_i4 + cosv(:) * fcut_i2 * fcut_i3 * dfcut_i4 * tmp_dxp(ix, i4+1)/r_central(i4)
                  end do
              end if
            else
              ! ---- multispecies modes 1/2/3: 4-body cluster (sa, sj1, sj2, sj3) channel routing ----
              sj1 = i_type_db(i2 + 1)
              sj2 = i_type_db(i3 + 1)
              sj3 = i_type_db(i4 + 1)
              call ftnbody_chem_route_4b(ftnbody_chem_mode, sa, sj1, sj2, sj3, nsp, nrff, nchan, ic, off, sgn, chem_w)
              gE(1:nrff) = cosv(1:nrff) * fcut_i2 * fcut_i3 * fcut_i4
              if (ftnbody_chem_mode == 2) then
                do irk = 1, nchan
                  off = (irk - 1)*nrff
                  tmp_real(off+1:off+nrff) = tmp_real(off+1:off+nrff) + chem_w(irk)*gE(1:nrff)
                end do
              else
                tmp_real(off+1:off+nrff) = tmp_real(off+1:off+nrff) + sgn*gE(1:nrff)
              end if
              if (desc_forces_local) then
                do ix = 1, 3
                  d2t(:) = -sinv(:)*Wmat(:,ix)
                  d3t(:) = -sinv(:)*Wmat(:,3+ix)
                  d4t(:) = -sinv(:)*Wmat(:,6+ix)
                  gF2(1:nrff) = d2t(1:nrff)*fcut_i2*fcut_i3*fcut_i4 + cosv(1:nrff)*dfcut_i2*fcut_i3*fcut_i4*tmp_dxp(ix, i2+1)/r_central(i2)
                  gF3(1:nrff) = d3t(1:nrff)*fcut_i2*fcut_i3*fcut_i4 + cosv(1:nrff)*fcut_i2*dfcut_i3*fcut_i4*tmp_dxp(ix, i3+1)/r_central(i3)
                  gF4(1:nrff) = d4t(1:nrff)*fcut_i2*fcut_i3*fcut_i4 + cosv(1:nrff)*fcut_i2*fcut_i3*dfcut_i4*tmp_dxp(ix, i4+1)/r_central(i4)
                  if (ftnbody_chem_mode == 2) then
                    do irk = 1, nchan
                      off = (irk - 1)*nrff
                      tmp_dreal(off+1:off+nrff,i2,ix) = tmp_dreal(off+1:off+nrff,i2,ix) + chem_w(irk)*gF2(1:nrff)
                      tmp_dreal(off+1:off+nrff,i3,ix) = tmp_dreal(off+1:off+nrff,i3,ix) + chem_w(irk)*gF3(1:nrff)
                      tmp_dreal(off+1:off+nrff,i4,ix) = tmp_dreal(off+1:off+nrff,i4,ix) + chem_w(irk)*gF4(1:nrff)
                    end do
                  else
                    tmp_dreal(off+1:off+nrff,i2,ix) = tmp_dreal(off+1:off+nrff,i2,ix) + sgn*gF2(1:nrff)
                    tmp_dreal(off+1:off+nrff,i3,ix) = tmp_dreal(off+1:off+nrff,i3,ix) + sgn*gF3(1:nrff)
                    tmp_dreal(off+1:off+nrff,i4,ix) = tmp_dreal(off+1:off+nrff,i4,ix) + sgn*gF4(1:nrff)
                  end if
                end do
              end if
            end if
  
            if (debug_time) then 
              time33 = MY_MPI_WTIME()
              t4b_inner_init =  t4b_inner_init + time11 - time00 
              t4b_inner_desc = t4b_inner_desc + time22 - time11 
              t4b_inner_deriv =  t4b_inner_deriv + time33 - time22 
            end if   

    
          end do ! end_i4 
        end do   ! end_i3
      end do     ! end_i2 
      !norm = norm/2 

      tmp_real(:) = tmp_real(:)*norm 

      if (debug_time) time22 = MY_MPI_WTIME()

      !cos! Nos this part is obsolote ... 
      !cos!   if (desc_forces_local) then 
      !cos! 
      !cos!   do kk = 1, max_neigh_local  
      !cos!     rtmp(:,:) = 0.d0 
      !cos! 
      !cos!     do i3 = kk+1,max_neigh_local
      !cos!       do i4 = i3 +1, max_neigh_local
      !cos!         ii = innvec(kk,i3,i4)
      !cos!         do ix = 1,3
      !cos!         rtmp(:,ix) = rtmp(:,ix) + d2_tmpr(:, ix, ii)
      !cos!         end do 
      !cos!       end do 
      !cos!     end do     
      !cos! 
      !cos!     do i2 = 1, max_neigh_local
      !cos!       if (kk > i2) then 
      !cos!       do i4 = kk+1, max_neigh_local
      !cos!         ii = innvec(i2,kk,i4)
      !cos!         do ix =1,3
      !cos!         !rtmp(:,ix) = rtmp(:,ix) + d3_tmpr(i2,kk,i4)
      !cos!         rtmp(:,ix) = rtmp(:,ix) + d3_tmpr(:, ix, ii)
      !cos!         end do 
      !cos!       end do 
      !cos!       end if 
      !cos!     end do     
      !cos! 
      !cos!     do i2 = 1, max_neigh_local
      !cos!       if (kk > (i2 +1) ) then 
      !cos!       do i3 = i2+1, max_neigh_local
      !cos!         if (kk > i3 ) then 
      !cos!         ii = innvec(i2,i3,kk)
      !cos!         do ix =1,3 
      !cos!         !rtmp(:,ix) = rtmp(:,ix) + d4_tmpr(i2,i3,kk)
      !cos!         rtmp(:,ix) = rtmp(:,ix) + d4_tmpr(:, ix, ii)
      !cos!         end do 
      !cos!         end if 
      !cos!       end do
      !cos!       end if  
      !cos!     end do  
      !cos! 
      !cos!     do ix =1, 3
      !cos!       tmp_dreal(:,kk,ix) = rtmp(:,ix)*norm  
      !cos!     end do       
      !cos!   end do   ! end kk 
      !cos! 
      !cos! end if
      if (desc_forces_local) then 
        do kk = 1, max_neigh_local
          do ix =1, 3
            tmp_dreal(:,kk,ix) = tmp_dreal(:,kk,ix)*norm  
          end do       
        end do   ! end kk 
      end if  
      
      if (debug_time) then 
        time33 = MY_MPI_WTIME()
        t4b_inner_last =  t4b_inner_last + time33 - time22 
      end if 

      deallocate(tmpf, chem_w)
     end subroutine ftnbody_order_4

    ! ------------------------------------------------------------------
    ! Channel routing of one 4-body cluster (central species sa, neighbour
    ! species sj1, sj2, sj3) for each ftnbody_chem_mode.  The compact 4-body q is
    ! permutation invariant in the 3 neighbours, so exact/hash channels use the
    ! sorted neighbour multiset and the low-rank factor is symmetrized over the 6
    ! neighbour permutations (S_3 permanent).  See docs/perspective_ftnbody.md.
    ! ------------------------------------------------------------------
    subroutine ftnbody_chem_route_4b(mode, sa, sj1, sj2, sj3, nsp, nrff, nchan, ic, off, sgn, chem_w)
      use module_ftnbody, only: ftnbody_chem_center, ftnbody_chem_neigh
      implicit none
      integer, intent(in)  :: mode, sa, sj1, sj2, sj3, nsp, nrff, nchan
      integer, intent(out) :: ic, off
      real(kind_double), intent(out) :: sgn
      real(kind_double), dimension(:), intent(inout) :: chem_w
      integer :: irk, p, q, t, a, b, mtri, mset
      integer(8) :: hkey
      real(kind_double) :: e1a, e1b, e1c, e2a, e2b, e2c, e3a, e3b, e3c
      sgn = 1.d0
      ic  = 1
      ! sort the neighbour triple p <= q <= t
      p = sj1; q = sj2; t = sj3
      if (p > q) call iswap(p, q)
      if (q > t) call iswap(q, t)
      if (p > q) call iswap(p, q)
      select case (mode)
      case (1)                 ! exact: central x unordered neighbour multiset of size 3
        mtri = nsp*(nsp+1)*(nsp+2)/6                       ! # size-3 multisets
        mset = 1
        do a = 1, p - 1
          mset = mset + (nsp - a + 1)*(nsp - a + 2)/2      ! size-2 multisets from {a..S}
        end do
        do b = p, q - 1
          mset = mset + (nsp - b + 1)                      ! size-1 multisets from {b..S}
        end do
        mset = mset + (t - q)
        ic = (sa - 1)*mtri + mset
      case (2)                 ! low-rank, symmetrized over the 6 neighbour permutations (S_3 permanent)
        do irk = 1, nchan
          e1a = ftnbody_chem_neigh(sj1, irk, 1, 4); e2a = ftnbody_chem_neigh(sj1, irk, 2, 4); e3a = ftnbody_chem_neigh(sj1, irk, 3, 4)
          e1b = ftnbody_chem_neigh(sj2, irk, 1, 4); e2b = ftnbody_chem_neigh(sj2, irk, 2, 4); e3b = ftnbody_chem_neigh(sj2, irk, 3, 4)
          e1c = ftnbody_chem_neigh(sj3, irk, 1, 4); e2c = ftnbody_chem_neigh(sj3, irk, 2, 4); e3c = ftnbody_chem_neigh(sj3, irk, 3, 4)
          chem_w(irk) = ftnbody_chem_center(sa, irk, 4) * &
               ( e1a*e2b*e3c + e1a*e2c*e3b + e1b*e2a*e3c &
               + e1b*e2c*e3a + e1c*e2a*e3b + e1c*e2b*e3a )
        end do
        ic = 1
      case (3)                 ! hash of (sa, sorted triple)
        hkey = ((int(sa-1,8)*int(nsp,8) + int(p-1,8))*int(nsp,8) + int(q-1,8))*int(nsp,8) + int(t-1,8)
        ic  = int(mod(hkey*2654435761_8, int(nchan, 8)), kind=4) + 1
        if (mod(hkey*40503_8 + 12345_8, 2_8) == 0_8) then
          sgn = 1.d0
        else
          sgn = -1.d0
        end if
      end select
      off = (ic - 1)*nrff
    contains
      subroutine iswap(m, n)
        integer, intent(inout) :: m, n
        integer :: k
        k = m; m = n; n = k
      end subroutine iswap
    end subroutine ftnbody_chem_route_4b

     subroutine ftnbody_order_5(type_db_ja, max_neigh_local, i_type_db, i_central, r_central, tmp_dxp, ur_central, d_ur_central, desc_forces_local)
     
      ! compute all the body4 triangle with the following structure. Please pay attention that the distances
      ! are already r_cut -ed.
      !
      !     i4
      !    /
      !  x3
      ! / 
      ! j--x1--i2
      !  
      ! | x2
      ! x4  \
      ! |    i3
      ! i5
      !
      ! Angles are defined like so:
      ! x5: i2-i3
      ! x6: i2-i4
      ! x7: i3-i4
      ! x8: i2-i5
      ! x9: i3-i5
      ! x10: i4-i5
      !
       use module_ftnbody, only: tmp_real, tmp_dreal, ftnbody_omega_5b, ftnbody_phase_random_5b, delta_rff, dim_rff, &
                                 init_mode_ftnbody, covar_matrix_5b , count_5b, dim_qbody, mean_5b, ltmp_mean_ftnbody, &
                                 r_cut_ft5b, r_cut_width_ft5b, &
                                 ftnbody_chem_mode, ftnbody_n_channels, ftnbody_chem_center, ftnbody_chem_neigh
       use module_chemical_species, only: fix_no_of_elements
       use time_check_general, only : debug_time, MY_MPI_WTIME
       use module_body_desc, only : dim_desc_body
       use module_neigh_local, only: r_cut_in, r_cut_width_in, r_cut_pair_in, fcut_rij_inout
       implicit none
       integer, intent(in) ::  type_db_ja,  max_neigh_local
       integer, dimension(:), intent(in)  :: i_type_db, i_central
       real(kind_double), dimension(:), intent(in) :: r_central, ur_central 
       real(kind_double), dimension(:,:), intent(in) :: tmp_dxp
       real(kind_double), dimension(:,:), intent(in) :: d_ur_central
       logical, intent(in) :: desc_forces_local 
       integer :: i2, i3, i4, i5, ix, kk, ii,  dimf
       integer :: iiii, jjjj
       real(kind_double), dimension(:), allocatable :: tmpf
       real(kind_double), dimension(:,:), allocatable :: ttmp 
       real(kind_double) :: norm, one=1.d0 ,  time00, time11, time22
       integer :: incx = 1, incy = 1, icount_nn
       real(kind_double), dimension(3)   :: d2_x1, d2_x2, d2_x3, d2_x4, d2_x5, d2_x6, d2_x7, d2_x8, d2_x9, d2_x10, &
                                            d3_x1, d3_x2, d3_x3, d3_x4, d3_x5, d3_x6, d3_x7, d3_x8, d3_x9, d3_x10, &
                                            d4_x1, d4_x2, d4_x3, d4_x4, d4_x5, d4_x6, d4_x7, d4_x8, d4_x9, d4_x10, &
                                            d5_x1, d5_x2, d5_x3, d5_x4, d5_x5, d5_x6, d5_x7, d5_x8, d5_x9, d5_x10
       real(kind_double)        :: cos_2j3, cos_2j4, cos_2j5, cos_3j4, cos_3j5, cos_4j5, &
                                   x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, &
                                   x1_2, x1_3, x1_4, &
                                  x2_2, x2_3, x2_4, &
                                  x3_2, x3_3, x3_4, &
                                  x4_2, x4_3, x4_4, &
                                  x5_2, x5_3, x5_4, x5_5, &
                                  x6_2, x6_3, x6_4, x6_5, &
                                  x7_2, x7_3, x7_4, x7_5, &
                                  x8_2, x8_3, x8_4, x8_5, &
                                  x9_2, x9_3, x9_4, x9_5, &
                                  x10_2, x10_3, x10_4, x10_5
       integer, parameter :: dimq = 10
       real(kind_double), dimension(dimq)   :: qsym
       real(kind_double), dimension(3, dimq)      :: d2_qsym, d3_qsym, d4_qsym, d5_qsym
       real(kind_double), dimension(dim_rff(5)) :: d2t, d3t, d4t, d5t
       ! hoisted cos/sin of the Fourier phase (evaluated once per cluster)
       real(kind_double), dimension(dim_rff(5)) :: cosv, sinv
       ! R3: slot contraction as one dgemm  Wmat = omega^T * DQmat,
       ! DQmat(i, (t-1)*3+ix) = d{t}_qsym(ix, i); -sin applied afterwards
       real(kind_double), dimension(dimq, 12) :: DQmat
       real(kind_double), dimension(dim_rff(5), 12) :: Wmat
       real(kind_double) :: zero = 0.d0
       ! --- multispecies chemistry (docs/perspective_ftnbody.md) ---
       integer :: nrff, nsp, nchan, sa, sj1, sj2, sj3, sj4, ic, off, irk
       real(kind_double) :: sgn
       real(kind_double), dimension(dim_rff(5)) :: gE, gF2, gF3, gF4, gF5
       real(kind_double), dimension(:), allocatable :: chem_w
       real(kind_double), dimension(:), allocatable :: fcut_all, dfcut_all
       real(kind_double) :: fcut, dfcut, fcut_i2, dfcut_i2, fcut_i3, dfcut_i3, fcut_i4, dfcut_i4, fcut_i5, dfcut_i5                               
       real(kind_double) :: local_r_cut_in
       integer :: type_fcut_in, type_fcut_out
       integer :: izozo 
       integer, dimension(size(i_type_db)) :: iv1_zozo  
       integer, dimension(size(i_central)) :: iv2_zozo  
 
       izozo =  type_db_ja
       iv1_zozo = i_type_db 
       iv2_zozo = i_central
       !TODOftnbody 
       type_fcut_in = 3 
       type_fcut_out = 2 
     
       nrff  = dim_rff(5)
       nsp   = fix_no_of_elements
       nchan = ftnbody_n_channels(5)
       sa    = type_db_ja

       allocate(tmpf(nrff))
       allocate(chem_w(max(nchan,1)))
       allocate(ttmp(nrff,dimq))
       allocate(fcut_all(max_neigh_local))
       allocate(dfcut_all(max_neigh_local))

       ! normalization uses the geometry feature count (dim_rff), not the expanded block
       norm = sqrt(1.d0)/sqrt(dble(nrff))*delta_rff(5)
       ! compute_qja_a_partir_de_rcentral
       dimf = nrff
       tmpf(:) = ftnbody_phase_random_5b
 
       if (max_neigh_local <= 1) return
       tmp_real(:) = 0.d0
       if (desc_forces_local) then
          ! dim = 1,3
          d2_x1(:) = 0.d0
          d2_x2(:) = 0.d0
          d2_x3(:) = 0.d0
          d2_x4(:) = 0.d0
          d2_x5(:) = 0.d0
          d2_x6(:) = 0.d0
          d2_x7(:) = 0.d0
          d2_x8(:) = 0.d0
          d2_x9(:) = 0.d0
          d2_x10(:) = 0.d0
     
          d3_x1(:) = 0.d0
          d3_x2(:) = 0.d0
          d3_x3(:) = 0.d0
          d3_x4(:) = 0.d0
          d3_x5(:) = 0.d0
          d3_x6(:) = 0.d0
          d3_x7(:) = 0.d0
          d3_x8(:) = 0.d0
          d3_x9(:) = 0.d0
          d3_x10(:) = 0.d0

          d4_x1(:) = 0.d0
          d4_x2(:) = 0.d0
          d4_x3(:) = 0.d0
          d4_x4(:) = 0.d0
          d4_x5(:) = 0.d0
          d4_x6(:) = 0.d0
          d4_x7(:) = 0.d0
          d4_x8(:) = 0.d0
          d4_x9(:) = 0.d0
          d4_x10(:) = 0.d0
          
          d5_x1(:) = 0.d0
          d5_x2(:) = 0.d0
          d5_x3(:) = 0.d0
          d5_x4(:) = 0.d0
          d5_x5(:) = 0.d0
          d5_x6(:) = 0.d0
          d5_x7(:) = 0.d0
          d5_x8(:) = 0.d0
          d5_x9(:) = 0.d0
          d5_x10(:) = 0.d0
          tmp_dreal(:, :, :) = 0.d0
       end if
 
       do ii = 1, max_neigh_local
          ! Pair-specific inner cutoff: i_type_db(ii+1) due to assumed-shape remapping (index 1 = central)
          local_r_cut_in = r_cut_pair_in(type_db_ja, i_type_db(ii + 1))
          call fcut_rij_inout(type_fcut_in, type_fcut_out, r_central(ii), r_cut_ft5b, r_cut_width_ft5b, local_r_cut_in, r_cut_width_in, desc_forces_local, fcut, dfcut)
          fcut_all(ii) = fcut
          dfcut_all(ii) = dfcut
       end do

       icount_nn = 0 
       do i2 = 1, max_neigh_local
         fcut_i2 = fcut_all(i2)
         if (fcut_i2 == 0.d0) cycle
         dfcut_i2 = dfcut_all(i2)

         x1 = ur_central(i2)
         x1_2 = x1**2
         x1_3 = x1_2*x1
         x1_4 = x1_2**2
         if (desc_forces_local) then
           d2_x1(1:3) = d_ur_central(1:3, i2)
         end if
         do i3 = i2 + 1, max_neigh_local
           fcut_i3 = fcut_all(i3)
           if (fcut_i3 == 0.d0) cycle
           dfcut_i3 = dfcut_all(i3)
           
           x2 = ur_central(i3)
           x2_2 = x2**2
           x2_3 = x2_2*x2
           x2_4 = x2_2**2

           cos_2j3 = dot_product(tmp_dxp(:, i2+1), tmp_dxp(:, i3+1))/(r_central(i2)*r_central(i3))
           x5 = cos_2j3
           x5_2 = x5**2
           x5_3 = x5_2*x5
           x5_4 = x5_2**2
           x5_5 = x5_4*x5
           if (desc_forces_local) then
             d3_x2(1:3) = d_ur_central(1:3, i3)
             d2_x5(1:3) = tmp_dxp(1:3, i3+1)/(r_central(i2)*r_central(i3)) - tmp_dxp(1:3, i2+1)*cos_2j3/r_central(i2)**2
             d3_x5(1:3) = tmp_dxp(1:3, i2+1)/(r_central(i2)*r_central(i3)) - tmp_dxp(1:3, i3+1)*cos_2j3/r_central(i3)**2
           end if
     
           ttmp(:,:)= 0.d0 
           do i4 = i3 + 1, max_neigh_local
              fcut_i4 = fcut_all(i4)
              if (fcut_i4 == 0.d0) cycle
              dfcut_i4 = dfcut_all(i4)
              
              x3 = ur_central(i4)
              x3_2 = x3**2
              x3_3 = x3_2*x3
              x3_4 = x3_2**2
            
              cos_2j4 = dot_product(tmp_dxp(:, i2+1), tmp_dxp(:, i4+1))/(r_central(i2)*r_central(i4))
              x6 = cos_2j4
              x6_2 = x6**2
              x6_3 = x6_2*x6
              x6_4 = x6_2**2
              x6_5 = x6_4*x6

              cos_3j4 = dot_product(tmp_dxp(:, i3+1), tmp_dxp(:, i4+1))/(r_central(i3)*r_central(i4))
              x7 = cos_3j4    
              x7_2 = x7**2
              x7_3 = x7_2*x7
              x7_4 = x7_2**2
              x7_5 = x7_4*x7
             
              if (desc_forces_local) then
                d4_x3(1:3) = d_ur_central(1:3, i4)

                d2_x6(1:3) = tmp_dxp(1:3, i4+1)/(r_central(i2)*r_central(i4)) - tmp_dxp(1:3, i2+1)*cos_2j4/r_central(i2)**2
                d4_x6(1:3) = tmp_dxp(1:3, i2+1)/(r_central(i2)*r_central(i4)) - tmp_dxp(1:3, i4+1)*cos_2j4/r_central(i4)**2
                d3_x7(1:3) = tmp_dxp(1:3, i4+1)/(r_central(i3)*r_central(i4)) - tmp_dxp(1:3, i3+1)*cos_3j4/r_central(i3)**2
                d4_x7(1:3) = tmp_dxp(1:3, i3+1)/(r_central(i3)*r_central(i4)) - tmp_dxp(1:3, i4+1)*cos_3j4/r_central(i4)**2
              end if
              
              do i5 = i4 + 1, max_neigh_local
                fcut_i5 = fcut_all(i5)
                if (fcut_i5 == 0.d0) cycle
                dfcut_i5 = dfcut_all(i5)

                if (debug_time) time00 = MY_MPI_WTIME()
                x4 = ur_central(i5)
                x4_2 = x4**2
                x4_3 = x4_2*x4
                x4_4 = x4_2**2
                
                cos_2j5 = dot_product(tmp_dxp(:, i2+1), tmp_dxp(:, i5+1))/(r_central(i2)*r_central(i5))
                x8 = cos_2j5
                x8_2 = x8**2
                x8_3 = x8_2*x8
                x8_4 = x8_2**2
                x8_5 = x8_4*x8

                cos_3j5 = dot_product(tmp_dxp(:, i3+1), tmp_dxp(:, i5+1))/(r_central(i3)*r_central(i5))
                x9 = cos_3j5
                x9_2 = x9**2
                x9_3 = x9_2*x9
                x9_4 = x9_2**2
                x9_5 = x9_4*x9

                cos_4j5 = dot_product(tmp_dxp(:, i4+1), tmp_dxp(:, i5+1))/(r_central(i4)*r_central(i5))
                x10 = cos_4j5
                x10_2 = x10**2
                x10_3 = x10_2*x10
                x10_4 = x10_2**2
                x10_5 = x10_4*x10

                if (desc_forces_local) then
                  d5_x4(1:3)  = d_ur_central(1:3, i5)
                  d2_x8(1:3)  = tmp_dxp(1:3, i5+1)/(r_central(i2)*r_central(i5)) - tmp_dxp(1:3, i2+1)*cos_2j5/r_central(i2)**2
                  d5_x8(1:3)  = tmp_dxp(1:3, i2+1)/(r_central(i2)*r_central(i5)) - tmp_dxp(1:3, i5+1)*cos_2j5/r_central(i5)**2
                  d3_x9(1:3)  = tmp_dxp(1:3, i5+1)/(r_central(i3)*r_central(i5)) - tmp_dxp(1:3, i3+1)*cos_3j5/r_central(i3)**2
                  d5_x9(1:3)  = tmp_dxp(1:3, i3+1)/(r_central(i3)*r_central(i5)) - tmp_dxp(1:3, i5+1)*cos_3j5/r_central(i5)**2
                  d4_x10(1:3) = tmp_dxp(1:3, i5+1)/(r_central(i4)*r_central(i5)) - tmp_dxp(1:3, i4+1)*cos_4j5/r_central(i4)**2
                  d5_x10(1:3) = tmp_dxp(1:3, i4+1)/(r_central(i4)*r_central(i5)) - tmp_dxp(1:3, i5+1)*cos_4j5/r_central(i5)**2
                end if
                ! Polynomes are a bit different from which G.Csanyi paper used
                ! qsym(1) = x1 + x2 + x3 + x4
                ! qsym(2) = x4 + x5 + x6 + x7 + x8 + x9 + x10
                ! qsym(3) = x1_2 + x2_2 + x3_2 + x4_2
                ! qsym(4) = x4**2 + x5**2 + x6**2 + x7**2 + x8**2 + x9**2 + x10**2
                ! qsym(5) = x1*(x5+x6+x8) + x2*(x5+x7+x9) + x3*(x6+x7+x10) + x4*(x8+x9+x10)
                ! qsym(6) = x5*x6*x8 + x5*x7*x9 + x6*x7*x10 + x8*x9*x10
                ! qsym(7) = x1_3 + x2_3 + x3_3 + x4_3
                ! qsym(8) = x4**3 + x5**3 + x6**3 + x7**3 + x8**3 + x9**3 + x10**3 
                ! qsym(9) = x1_2*(x5+x6+x8) + x2_2*(x5+x7+x9) + x3_2*(x6+x7+x10) + x4_2*(x8+x9+x10)
                ! qsym(10) = x1_4 + x2_4 + x3_4 + x4_4

                qsym(1) = x1 + x2 + x3 + x4 
                qsym(2) = x5 + x6 + x7 + x8 + x9 + x10 
                qsym(3) = x1_2 + x2_2 + x3_2 + x4_2 
                qsym(4) = x5_2 + x6_2 + x7_2 + x8_2 + x9_2 + x10_2 
                qsym(5) = x1_3 + x2_3 + x3_3 + x4_3 
                qsym(6) = x5_3 + x6_3 + x7_3 + x8_3 + x9_3 + x10_3  
                qsym(7) = x1_4 + x2_4 + x3_4 + x4_4 
                qsym(8) = x5_4 + x6_4 + x7_4 + x8_4 + x9_4 + x10_4 
                qsym(9) = x5_5 + x6_5 + x7_5 + x8_5 + x9_5 + x10_5
                qsym(10) = x5*x6*x7*x8*x9*x10 

                if (init_mode_ftnbody) then
                  count_5b = count_5b + 1
                  if (ltmp_mean_ftnbody) then 
                   mean_5b(:) = mean_5b(:) + qsym(:)
                  else 
                    do iiii = 1, dim_qbody(5)
                      do jjjj = 1, dim_qbody(5) 
                        covar_matrix_5b(iiii,jjjj) = covar_matrix_5b(iiii,jjjj) + (qsym(iiii)  - mean_5b(iiii)) &
                                                                                 *(qsym(jjjj) - mean_5b(jjjj))
                      end do
                    end do
                  end if 
                  cycle
                end if

                if (desc_forces_local) then
                  ! Look a sqym(5) or qsym(9) to know wich ii is in which angle
                  ! 2 is in 1 bond ; 5,6,8 - angles
                  d2_qsym(1:3, 1) = d2_x1(1:3)
                  d2_qsym(1:3, 2) = d2_x5(1:3) + d2_x6(1:3) + d2_x8(1:3)
                  d2_qsym(1:3, 3) = 2.d0*x1*d2_x1(1:3)
                  d2_qsym(1:3, 4) = 2.d0*(d2_x5(1:3)*x5 + d2_x6(1:3)*x6 + d2_x8(1:3)*x8)
                  d2_qsym(1:3, 5) = 3.d0*x1_2*d2_x1(1:3)
                  d2_qsym(1:3, 6) = 3.d0*(d2_x5(1:3)*x5_2 + d2_x6(1:3)*x6_2 + d2_x8(1:3)*x8_2) 
                  d2_qsym(1:3, 7) = 4.d0*x1_3*d2_x1(1:3)
                  d2_qsym(1:3, 8) = 4.d0*(d2_x5(1:3)*x5_3 + d2_x6(1:3)*x6_3 + d2_x8(1:3)*x8_3)
                  d2_qsym(1:3, 9) = 5.d0*(d2_x5(1:3)*x5_4 + d2_x6(1:3)*x6_4 + d2_x8(1:3)*x8_4)
                  d2_qsym(1:3, 10) = x7*x9*x10*(x6*x8*d2_x5(1:3) + x5*x8*d2_x6(1:3) + x5*x6*d2_x8(1:3))
      
                  ! 3 is in 2 bond ; 5,7,9 - angles
                  d3_qsym(1:3, 1) = d3_x2(1:3)
                  d3_qsym(1:3, 2) = d3_x5(1:3) + d3_x7(1:3) + d3_x9(1:3)
                  d3_qsym(1:3, 3) = 2.d0*x2*d3_x2(1:3)
                  d3_qsym(1:3, 4) = 2.d0*(d3_x5(1:3)*x5 + d3_x7(1:3)*x7 + d3_x9(1:3)*x9)
                  d3_qsym(1:3, 5) = 3.d0*x2_2*d3_x2(1:3)
                  d3_qsym(1:3, 6) = 3.d0*(d3_x5(1:3)*x5_2 + d3_x7(1:3)*x7_2 + d3_x9(1:3)*x9_2)
                  d3_qsym(1:3, 7) = 4.d0*x2_3*d3_x2(1:3)
                  d3_qsym(1:3, 8) = 4.d0*(d3_x5(1:3)*x5_3 + d3_x7(1:3)*x7_3 + d3_x9(1:3)*x9_3)
                  d3_qsym(1:3, 9) = 5.d0*(d3_x5(1:3)*x5_4 + d3_x7(1:3)*x7_4 + d3_x9(1:3)*x9_4)
                  d3_qsym(1:3, 10) = x6*x8*x10*(x7*x9*d3_x5(1:3) + x5*x9*d3_x7(1:3) + x5*x7*d3_x9(1:3))
                  
                  ! 4 is in 3 bond ; 6,7,10 - angles
                  d4_qsym(1:3, 1) = d4_x3(1:3)
                  d4_qsym(1:3, 2) = d4_x6(1:3) + d4_x7(1:3) + d4_x10(1:3)
                  d4_qsym(1:3, 3) = 2.d0*x3*d4_x3(1:3)
                  d4_qsym(1:3, 4) = 2.d0*(d4_x6(1:3)*x6 + d4_x7(1:3)*x7 + d4_x10(1:3)*x10)
                  d4_qsym(1:3, 5) = 3.d0*x3_2*d4_x3(1:3)
                  d4_qsym(1:3, 6) = 3.d0*(d4_x6(1:3)*x6_2 + d4_x7(1:3)*x7_2 + d4_x10(1:3)*x10_2)
                  d4_qsym(1:3, 7) = 4.d0*x3_3*d4_x3(1:3)
                  d4_qsym(1:3, 8) = 4.d0*(d4_x6(1:3)*x6_3 + d4_x7(1:3)*x7_3 + d4_x10(1:3)*x10_3)
                  d4_qsym(1:3, 9) = 5.d0*(d4_x6(1:3)*x6_4 + d4_x7(1:3)*x7_4 + d4_x10(1:3)*x10_4)
                  d4_qsym(1:3, 10) = x5*x8*x9*(x6*x7*d4_x10(1:3) + x6*x10*d4_x7(1:3) + x7*x10*d4_x6(1:3))

                  ! 5 is in 4 bond ; 8,9,10 - angles
                  d5_qsym(1:3, 1) = d5_x4(1:3)
                  d5_qsym(1:3, 2) = d5_x8(1:3) + d5_x9(1:3) + d5_x10(1:3)
                  d5_qsym(1:3, 3) = 2.d0*x4*d5_x4(1:3)
                  d5_qsym(1:3, 4) = 2.d0*(d5_x8(1:3)*x8 + d5_x9(1:3)*x9 + d5_x10(1:3)*x10)
                  d5_qsym(1:3, 5) = 3.d0*x4_2*d5_x4
                  d5_qsym(1:3, 6) = 3.d0*(d5_x8(1:3)*x8_2 + d5_x9(1:3)*x9_2 + d5_x10(1:3)*x10_2)
                  d5_qsym(1:3, 7) = 4.d0*x4_3*d5_x4(1:3)
                  d5_qsym(1:3, 8) = 4.d0*(d5_x8(1:3)*x8_3 + d5_x9(1:3)*x9_3 + d5_x10(1:3)*x10_3)
                  d5_qsym(1:3, 9) = 5.d0*(d5_x8(1:3)*x8_4 + d5_x9(1:3)*x9_4 + d5_x10(1:3)*x10_4)
                  d5_qsym(1:3, 10) = x5*x6*x7*(x8*x9*d5_x10(1:3) + x8*x10*d5_x9(1:3) + x9*x10*d5_x8(1:3))
                                    !I DEBUG UNTIL HERE ......................!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
              end if
                if (debug_time) time11 = MY_MPI_WTIME()
  
                tmpf(:) = ftnbody_phase_random_5b(:)
                call dgemv('T', dimq, dimf, one, ftnbody_omega_5b, dimq, qsym, incx, one, tmpf, incy)

                if (debug_time) time22 = MY_MPI_WTIME()

                ! evaluate the transcendentals once per cluster; R3: the slot
                ! contraction is one dgemm on the raw omega, -sin applied after
                cosv(:) = cos(tmpf(:))
                if (desc_forces_local) then
                  sinv(:) = sin(tmpf(:))
                  do ii = 1, dimq
                    DQmat(ii, 1:3)   = d2_qsym(1:3, ii)
                    DQmat(ii, 4:6)   = d3_qsym(1:3, ii)
                    DQmat(ii, 7:9)   = d4_qsym(1:3, ii)
                    DQmat(ii, 10:12) = d5_qsym(1:3, ii)
                  end do
                  call dgemm('T', 'N', nrff, 12, dimq, one, ftnbody_omega_5b, dimq, DQmat, dimq, zero, Wmat, nrff)
                end if

                if (ftnbody_chem_mode == 0) then
                ! ---- single species: bit-identical to the original implementation ----
                tmp_real(:) = tmp_real(:) + cosv(:) * fcut_i2 * fcut_i3 * fcut_i4 * fcut_i5
                if (desc_forces_local) then
                  do ix = 1, 3
                    d2t(:) = -sinv(:)*Wmat(:,ix)
                    d3t(:) = -sinv(:)*Wmat(:,3+ix)
                    d4t(:) = -sinv(:)*Wmat(:,6+ix)
                    d5t(:) = -sinv(:)*Wmat(:,9+ix)

                    tmp_dreal(:,i2,ix) = tmp_dreal(:,i2,ix) + d2t(:) * fcut_i2 * fcut_i3 * fcut_i4 * fcut_i5 &
                                          +  cosv(:) * dfcut_i2 * fcut_i3 * fcut_i4 * fcut_i5 * tmp_dxp(ix, i2+1)/r_central(i2)
                    tmp_dreal(:,i3,ix) = tmp_dreal(:,i3,ix) + d3t(:) * fcut_i2 * fcut_i3 * fcut_i4 * fcut_i5 &
                                          +  cosv(:) * fcut_i2 * dfcut_i3 * fcut_i4 * fcut_i5 * tmp_dxp(ix, i3+1)/r_central(i3)
                    tmp_dreal(:,i4,ix) = tmp_dreal(:,i4,ix) + d4t(:) * fcut_i2 * fcut_i3 * fcut_i4 * fcut_i5 &
                                          +  cosv(:) * fcut_i2 * fcut_i3 * dfcut_i4 * fcut_i5 * tmp_dxp(ix, i4+1)/r_central(i4)
                    tmp_dreal(:,i5,ix) = tmp_dreal(:,i5,ix) + d5t(:) * fcut_i2 * fcut_i3 * fcut_i4 * fcut_i5 &
                                          +  cosv(:) * fcut_i2 * fcut_i3 * fcut_i4 * dfcut_i5 * tmp_dxp(ix, i5+1)/r_central(i5)
                    end do
                end if
                else
                ! ---- multispecies modes 1/2/3: 5-body cluster (sa, sj1..sj4) channel routing ----
                sj1 = i_type_db(i2 + 1)
                sj2 = i_type_db(i3 + 1)
                sj3 = i_type_db(i4 + 1)
                sj4 = i_type_db(i5 + 1)
                call ftnbody_chem_route_5b(ftnbody_chem_mode, sa, sj1, sj2, sj3, sj4, nsp, nrff, nchan, ic, off, sgn, chem_w)
                gE(1:nrff) = cosv(1:nrff) * fcut_i2 * fcut_i3 * fcut_i4 * fcut_i5
                if (ftnbody_chem_mode == 2) then
                  do irk = 1, nchan
                    off = (irk - 1)*nrff
                    tmp_real(off+1:off+nrff) = tmp_real(off+1:off+nrff) + chem_w(irk)*gE(1:nrff)
                  end do
                else
                  tmp_real(off+1:off+nrff) = tmp_real(off+1:off+nrff) + sgn*gE(1:nrff)
                end if
                if (desc_forces_local) then
                  do ix = 1, 3
                    d2t(:) = -sinv(:)*Wmat(:,ix)
                    d3t(:) = -sinv(:)*Wmat(:,3+ix)
                    d4t(:) = -sinv(:)*Wmat(:,6+ix)
                    d5t(:) = -sinv(:)*Wmat(:,9+ix)
                    gF2(1:nrff) = d2t(1:nrff)*fcut_i2*fcut_i3*fcut_i4*fcut_i5 + cosv(1:nrff)*dfcut_i2*fcut_i3*fcut_i4*fcut_i5*tmp_dxp(ix, i2+1)/r_central(i2)
                    gF3(1:nrff) = d3t(1:nrff)*fcut_i2*fcut_i3*fcut_i4*fcut_i5 + cosv(1:nrff)*fcut_i2*dfcut_i3*fcut_i4*fcut_i5*tmp_dxp(ix, i3+1)/r_central(i3)
                    gF4(1:nrff) = d4t(1:nrff)*fcut_i2*fcut_i3*fcut_i4*fcut_i5 + cosv(1:nrff)*fcut_i2*fcut_i3*dfcut_i4*fcut_i5*tmp_dxp(ix, i4+1)/r_central(i4)
                    gF5(1:nrff) = d5t(1:nrff)*fcut_i2*fcut_i3*fcut_i4*fcut_i5 + cosv(1:nrff)*fcut_i2*fcut_i3*fcut_i4*dfcut_i5*tmp_dxp(ix, i5+1)/r_central(i5)
                    if (ftnbody_chem_mode == 2) then
                      do irk = 1, nchan
                        off = (irk - 1)*nrff
                        tmp_dreal(off+1:off+nrff,i2,ix) = tmp_dreal(off+1:off+nrff,i2,ix) + chem_w(irk)*gF2(1:nrff)
                        tmp_dreal(off+1:off+nrff,i3,ix) = tmp_dreal(off+1:off+nrff,i3,ix) + chem_w(irk)*gF3(1:nrff)
                        tmp_dreal(off+1:off+nrff,i4,ix) = tmp_dreal(off+1:off+nrff,i4,ix) + chem_w(irk)*gF4(1:nrff)
                        tmp_dreal(off+1:off+nrff,i5,ix) = tmp_dreal(off+1:off+nrff,i5,ix) + chem_w(irk)*gF5(1:nrff)
                      end do
                    else
                      tmp_dreal(off+1:off+nrff,i2,ix) = tmp_dreal(off+1:off+nrff,i2,ix) + sgn*gF2(1:nrff)
                      tmp_dreal(off+1:off+nrff,i3,ix) = tmp_dreal(off+1:off+nrff,i3,ix) + sgn*gF3(1:nrff)
                      tmp_dreal(off+1:off+nrff,i4,ix) = tmp_dreal(off+1:off+nrff,i4,ix) + sgn*gF4(1:nrff)
                      tmp_dreal(off+1:off+nrff,i5,ix) = tmp_dreal(off+1:off+nrff,i5,ix) + sgn*gF5(1:nrff)
                    end if
                    end do
                end if
                end if
    
              !if (debug_time) then 
                !time33 = MY_MPI_WTIME()
                !t5b_inner_init =  t5b_inner_init + time11 - time00
                !t5b_inner_desc = t5b_inner_desc + time22 - time11
                !t5b_inner_deriv =  t5b_inner_deriv + time33 - time22
              !end if   
  
            end do ! end_i5
           end do ! end_i4 
         end do   ! end_i3
       end do     ! end_i2 

       tmp_real(:) = tmp_real(:)*norm 
       if (desc_forces_local) then 
        do kk = 1, max_neigh_local
          do ix =1, 3
            tmp_dreal(:,kk,ix) = tmp_dreal(:,kk,ix)*norm  
          end do       
        end do   ! end kk 
      end if  



       deallocate(tmpf, chem_w)
      end subroutine ftnbody_order_5

    ! ------------------------------------------------------------------
    ! Channel routing of one 5-body cluster (central species sa, neighbour
    ! species sj1..sj4) for each ftnbody_chem_mode.  The compact 5-body q is
    ! permutation invariant in the 4 neighbours, so exact/hash channels use the
    ! sorted neighbour multiset and the low-rank factor is symmetrized over the
    ! 24 neighbour permutations (S_4 permanent, via Ryser).  See
    ! docs/perspective_ftnbody.md.
    ! ------------------------------------------------------------------
    subroutine ftnbody_chem_route_5b(mode, sa, sj1, sj2, sj3, sj4, nsp, nrff, nchan, ic, off, sgn, chem_w)
      use module_ftnbody, only: ftnbody_chem_center, ftnbody_chem_neigh
      implicit none
      integer, intent(in)  :: mode, sa, sj1, sj2, sj3, sj4, nsp, nrff, nchan
      integer, intent(out) :: ic, off
      real(kind_double), intent(out) :: sgn
      real(kind_double), dimension(:), intent(inout) :: chem_w
      integer :: irk, p, q, t, u, a, b, c, sv(4), mquad, mset
      integer(8) :: hkey
      real(kind_double) :: amat(4,4)
      sgn = 1.d0
      ic  = 1
      sv = (/ sj1, sj2, sj3, sj4 /)
      call isort4(sv)
      p = sv(1); q = sv(2); t = sv(3); u = sv(4)
      select case (mode)
      case (1)                 ! exact: central x unordered neighbour multiset of size 4
        mquad = nsp*(nsp+1)*(nsp+2)*(nsp+3)/24
        mset = 1
        do a = 1, p - 1
          mset = mset + (nsp-a+1)*(nsp-a+2)*(nsp-a+3)/6    ! size-3 multisets from {a..S}
        end do
        do b = p, q - 1
          mset = mset + (nsp-b+1)*(nsp-b+2)/2              ! size-2 multisets from {b..S}
        end do
        do c = q, t - 1
          mset = mset + (nsp-c+1)                          ! size-1 multisets from {c..S}
        end do
        mset = mset + (u - t)
        ic = (sa - 1)*mquad + mset
      case (2)                 ! low-rank, S_4 permanent of the role-embedding matrix
        do irk = 1, nchan
          amat(1,1)=ftnbody_chem_neigh(sj1,irk,1,5); amat(1,2)=ftnbody_chem_neigh(sj2,irk,1,5)
          amat(1,3)=ftnbody_chem_neigh(sj3,irk,1,5); amat(1,4)=ftnbody_chem_neigh(sj4,irk,1,5)
          amat(2,1)=ftnbody_chem_neigh(sj1,irk,2,5); amat(2,2)=ftnbody_chem_neigh(sj2,irk,2,5)
          amat(2,3)=ftnbody_chem_neigh(sj3,irk,2,5); amat(2,4)=ftnbody_chem_neigh(sj4,irk,2,5)
          amat(3,1)=ftnbody_chem_neigh(sj1,irk,3,5); amat(3,2)=ftnbody_chem_neigh(sj2,irk,3,5)
          amat(3,3)=ftnbody_chem_neigh(sj3,irk,3,5); amat(3,4)=ftnbody_chem_neigh(sj4,irk,3,5)
          amat(4,1)=ftnbody_chem_neigh(sj1,irk,4,5); amat(4,2)=ftnbody_chem_neigh(sj2,irk,4,5)
          amat(4,3)=ftnbody_chem_neigh(sj3,irk,4,5); amat(4,4)=ftnbody_chem_neigh(sj4,irk,4,5)
          chem_w(irk) = ftnbody_chem_center(sa, irk, 5) * perm4(amat)
        end do
        ic = 1
      case (3)                 ! hash of (sa, sorted quadruple)
        hkey = (((int(sa-1,8)*int(nsp,8) + int(p-1,8))*int(nsp,8) + int(q-1,8))*int(nsp,8) &
                + int(t-1,8))*int(nsp,8) + int(u-1,8)
        ic  = int(mod(hkey*2654435761_8, int(nchan, 8)), kind=4) + 1
        if (mod(hkey*40503_8 + 12345_8, 2_8) == 0_8) then
          sgn = 1.d0
        else
          sgn = -1.d0
        end if
      end select
      off = (ic - 1)*nrff
    contains
      subroutine isort4(v)
        integer, intent(inout) :: v(4)
        integer :: i, j, k
        do i = 1, 3
          do j = i+1, 4
            if (v(j) < v(i)) then
              k = v(i); v(i) = v(j); v(j) = k
            end if
          end do
        end do
      end subroutine isort4
      ! permanent of a 4x4 matrix via Ryser's formula (n=4 -> (-1)^n=+1)
      real(kind_double) function perm4(am) result(pm)
        real(kind_double), intent(in) :: am(4,4)
        real(kind_double) :: rowsum(4), prod
        integer :: s, i, jj, pc
        pm = 0.d0
        do s = 1, 15            ! nonempty subsets of {1,2,3,4}; empty set contributes 0
          pc = popcnt(s)
          do i = 1, 4
            rowsum(i) = 0.d0
            do jj = 1, 4
              if (btest(s, jj-1)) rowsum(i) = rowsum(i) + am(i, jj)
            end do
          end do
          prod = rowsum(1)*rowsum(2)*rowsum(3)*rowsum(4)
          if (mod(pc, 2) == 0) then     ! (-1)^|S|
            pm = pm + prod
          else
            pm = pm - prod
          end if
        end do
      end function perm4
    end subroutine ftnbody_chem_route_5b

    ! ==================================================================
    ! ftnbody_model = 'gramm'  (docs/MLT5 sec. "Permutation-invariant
    ! Gram-coordinate Fourier representation")
    !
    ! Ordered Gram coordinates q_G = (u_1..u_m, c_12, ..., c_{m-1,m}) with
    ! m = n-1 neighbours; permutation invariance is imposed AFTER the Fourier
    ! map by an exact average over the S_m orbit:
    !
    !   D_a^(n) = delta_rff/sqrt(Z_n) sum_clusters C_ja *
    !             (1/m!) sum_pi chi_pi cos( omega^T Pi_pi q_G + b )
    !
    ! chi_pi is 1 (mode 0), the exact/hashed channel weight (modes 1/3,
    ! permutation invariant -> routed once per cluster), or the jointly
    ! symmetrized low-rank chemistry product (mode 2, per permutation):
    !   chi_{r,pi} = E^c(sa,r) * prod_p E^p( s_{j_pi(p)}, r )
    ! ------------------------------------------------------------------
    ! Dispatch wrapper: selects the per-order module arrays.
    ! ==================================================================
    subroutine ftnbody_gram_order(nord, type_db_ja, max_neigh_local, i_type_db, r_central, tmp_dxp, &
                                  ur_central, d_ur_central, desc_forces_local)
      use module_ftnbody, only: ftnbody_omega_2b, ftnbody_omega_3b, ftnbody_omega_4b, ftnbody_omega_5b, &
                                ftnbody_phase_random_2b, ftnbody_phase_random_3b, &
                                ftnbody_phase_random_4b, ftnbody_phase_random_5b, &
                                r_cut_ft2b, r_cut_width_ft2b, r_cut_ft3b, r_cut_width_ft3b, &
                                r_cut_ft4b, r_cut_width_ft4b, r_cut_ft5b, r_cut_width_ft5b, &
                                count_2b, count_3b, count_4b, count_5b, &
                                mean_2b, mean_3b, mean_4b, mean_5b, &
                                covar_matrix_2b, covar_matrix_3b, covar_matrix_4b, covar_matrix_5b
      implicit none
      integer, intent(in) :: nord, type_db_ja, max_neigh_local
      integer, dimension(:), intent(in)  :: i_type_db
      real(kind_double), dimension(:), intent(in) :: r_central, ur_central
      real(kind_double), dimension(:,:), intent(in) :: tmp_dxp
      real(kind_double), dimension(:,:), intent(in) :: d_ur_central
      logical, intent(in) :: desc_forces_local

      select case (nord)
      case (2)
        call ftnbody_gram_generic(2, type_db_ja, max_neigh_local, i_type_db, r_central, tmp_dxp, &
                                  ur_central, d_ur_central, desc_forces_local, &
                                  ftnbody_omega_2b, ftnbody_phase_random_2b, &
                                  r_cut_ft2b, r_cut_width_ft2b, count_2b, mean_2b, covar_matrix_2b)
      case (3)
        call ftnbody_gram_generic(3, type_db_ja, max_neigh_local, i_type_db, r_central, tmp_dxp, &
                                  ur_central, d_ur_central, desc_forces_local, &
                                  ftnbody_omega_3b, ftnbody_phase_random_3b, &
                                  r_cut_ft3b, r_cut_width_ft3b, count_3b, mean_3b, covar_matrix_3b)
      case (4)
        call ftnbody_gram_generic(4, type_db_ja, max_neigh_local, i_type_db, r_central, tmp_dxp, &
                                  ur_central, d_ur_central, desc_forces_local, &
                                  ftnbody_omega_4b, ftnbody_phase_random_4b, &
                                  r_cut_ft4b, r_cut_width_ft4b, count_4b, mean_4b, covar_matrix_4b)
      case (5)
        call ftnbody_gram_generic(5, type_db_ja, max_neigh_local, i_type_db, r_central, tmp_dxp, &
                                  ur_central, d_ur_central, desc_forces_local, &
                                  ftnbody_omega_5b, ftnbody_phase_random_5b, &
                                  r_cut_ft5b, r_cut_width_ft5b, count_5b, mean_5b, covar_matrix_5b)
      end select
    end subroutine ftnbody_gram_order

    ! ------------------------------------------------------------------
    ! Order-generic Gram-coordinate worker (energy descriptor block,
    ! analytic force block, and mean/covariance accumulation for the
    ! Fourier-length estimation passes).
    ! ------------------------------------------------------------------
    subroutine ftnbody_gram_generic(nord, type_db_ja, max_neigh_local, i_type_db, r_central, tmp_dxp, &
                                    ur_central, d_ur_central, desc_forces_local, &
                                    omega, phase, rcut_n, rcutw_n, cnt_n, mean_n, covar_n)
      use module_ftnbody, only: tmp_real, tmp_dreal, delta_rff, dim_rff, &
                                init_mode_ftnbody, ltmp_mean_ftnbody, &
                                ftnbody_chem_mode, ftnbody_n_channels, &
                                ftnbody_chem_center, ftnbody_chem_neigh, &
                                gram_nperm, gram_perm, gram_qmap, gram_pair_pos
      use module_chemical_species, only: fix_no_of_elements
      use module_neigh_local, only: r_cut_in, r_cut_width_in, r_cut_pair_in, fcut_rij_inout
      implicit none
      integer, intent(in) :: nord, type_db_ja, max_neigh_local
      integer, dimension(:), intent(in)  :: i_type_db
      real(kind_double), dimension(:), intent(in) :: r_central, ur_central
      real(kind_double), dimension(:,:), intent(in) :: tmp_dxp
      real(kind_double), dimension(:,:), intent(in) :: d_ur_central
      logical, intent(in) :: desc_forces_local
      real(kind_double), dimension(:,:), intent(in) :: omega   ! (dimq, nrff)
      real(kind_double), dimension(:),   intent(in) :: phase   ! (nrff)
      real(kind_double), intent(in)    :: rcut_n, rcutw_n
      real(kind_double), intent(inout) :: cnt_n
      real(kind_double), dimension(:),   intent(inout) :: mean_n
      real(kind_double), dimension(:,:), intent(inout) :: covar_n

      integer, parameter :: MAXM = 4, MAXQ = 10, MAXPER = 24
      real(kind_double), parameter :: one = 1.d0, zero = 0.d0
      integer :: m, dimq, nrff, nper, nsp, nchan, sa
      integer :: idx(MAXM), jp(MAXM), sj(MAXM)
      real(kind_double) :: rp(MAXM), up(MAXM), fcp(MAXM), dfcp(MAXM), cex(MAXM)
      real(kind_double) :: qg(MAXQ)
      real(kind_double) :: dq(3, MAXQ, MAXM)
      real(kind_double) :: cprod, cpq, winv, norm, sgn, dcomp, wsc, wchem
      integer :: p, q, t, i, ix, k, ipm, irk, ic, off, ii, iiii, jjjj, ipos
      integer :: type_fcut_in, type_fcut_out, incx, incy
      integer :: g, ngrp, icol0
      integer :: grp_of(MAXPER), grp_key(MAXPER)
      real(kind_double) :: local_r_cut_in, fcut, dfcut
      real(kind_double), dimension(:), allocatable :: fcut_all, dfcut_all
      real(kind_double), dimension(:), allocatable :: gF, chem_w
      ! R3: slot contraction as one dgemm per group  Wg = SJg(:,:,g) * DQmat,
      ! DQmat(i, (t-1)*3+ix) = dq(ix, i, t)
      real(kind_double) :: DQmat(MAXQ, 3*MAXM)
      real(kind_double), dimension(:,:), allocatable :: Wg
      ! orbit batching: all permuted phases evaluated in one dgemv + one cos/sin;
      ! omegaP holds the row-permuted omega so that theta_pi = omegaP_pi^T q_G
      real(kind_double), dimension(:,:), allocatable :: omegaP        ! (dimq, nrff*nper)
      real(kind_double), dimension(:,:,:), allocatable :: omegaPT     ! (nrff, dimq, nper)
      real(kind_double), dimension(:,:), allocatable :: theta_all, cosA, sinA   ! (nrff, nper)
      ! per-group orbit sums (permutations sharing a species arrangement):
      !   CEg(f,g)   =  sum_{pi in g} cos(theta_pi)
      !   SJg(f,j,g) = -sum_{pi in g} sin(theta_pi) * omega'_pi(j)
      real(kind_double), dimension(:,:), allocatable :: CEg           ! (nrff, nper)
      real(kind_double), dimension(:,:,:), allocatable :: SJg         ! (nrff, dimq, nper)
      real(kind_double), dimension(:,:), allocatable :: chem_wg       ! (nchan, nper)
      logical :: okcut

      m    = nord - 1
      dimq = m + m*(m - 1)/2
      nrff = dim_rff(nord)
      nper = gram_nperm(nord)
      winv = 1.d0/dble(nper)
      nsp  = fix_no_of_elements
      nchan = ftnbody_n_channels(nord)
      sa   = type_db_ja
      norm = 1.d0/sqrt(dble(nrff))*delta_rff(nord)
      incx = 1
      incy = 1

      tmp_real(:) = 0.d0
      if (desc_forces_local) tmp_dreal(:, :, :) = 0.d0
      if (max_neigh_local < m) return

      allocate(fcut_all(max_neigh_local), dfcut_all(max_neigh_local))
      allocate(gF(nrff))
      allocate(Wg(nrff, 3*MAXM))
      allocate(chem_w(max(nchan, 1)))
      allocate(omegaP(dimq, nrff*nper))
      allocate(omegaPT(nrff, dimq, nper))
      allocate(theta_all(nrff, nper), cosA(nrff, nper), sinA(nrff, nper))
      allocate(CEg(nrff, nper), SJg(nrff, dimq, nper))
      allocate(chem_wg(max(nchan, 1), nper))

      ! build, once per central atom, the row-permuted omega copies:
      ! (Pi_pi q)_i = q(gram_qmap(i)) => theta_pi = sum_j omegaP_pi(j) q(j)
      ! with omegaP_pi(gram_qmap(i), :) = omega(i, :)
      if (.not. init_mode_ftnbody) then
        do ipm = 1, nper
          icol0 = (ipm - 1)*nrff
          do i = 1, dimq
            omegaP(gram_qmap(i, ipm, nord), icol0+1:icol0+nrff) = omega(i, 1:nrff)
            omegaPT(1:nrff, gram_qmap(i, ipm, nord), ipm) = omega(i, 1:nrff)
          end do
        end do
      end if

      !TODOftnbody (same convention as the poly model)
      type_fcut_in  = 3
      type_fcut_out = 2
      do ii = 1, max_neigh_local
        ! Pair-specific inner cutoff: i_type_db(ii+1) due to assumed-shape remapping (index 1 = central)
        local_r_cut_in = r_cut_pair_in(type_db_ja, i_type_db(ii + 1))
        call fcut_rij_inout(type_fcut_in, type_fcut_out, r_central(ii), rcut_n, rcutw_n, &
                            local_r_cut_in, r_cut_width_in, desc_forces_local, fcut, dfcut)
        fcut_all(ii)  = fcut
        dfcut_all(ii) = dfcut
      end do

      ! enumerate every unordered cluster j_1 < j_2 < ... < j_m exactly once
      do p = 1, m
        idx(p) = p
      end do

      odometer: do
        okcut = .true.
        do p = 1, m
          if (fcut_all(idx(p)) == 0.d0) then
            okcut = .false.
            exit
          end if
        end do

        if (okcut) then
          ! ---- ordered Gram coordinates of the cluster ----
          cprod = 1.d0
          do p = 1, m
            jp(p)  = idx(p)
            sj(p)  = i_type_db(jp(p) + 1)
            rp(p)  = r_central(jp(p))
            up(p)  = ur_central(jp(p))
            fcp(p) = fcut_all(jp(p))
            dfcp(p) = dfcut_all(jp(p))
            cprod  = cprod*fcp(p)
            qg(p)  = up(p)
          end do
          do p = 1, m - 1
            do q = p + 1, m
              cpq = dot_product(tmp_dxp(1:3, jp(p)+1), tmp_dxp(1:3, jp(q)+1))/(rp(p)*rp(q))
              qg(gram_pair_pos(m, p, q)) = cpq
            end do
          end do

          if (init_mode_ftnbody) then
            ! Fourier-length estimation passes (mean, then covariance)
            cnt_n = cnt_n + 1
            if (ltmp_mean_ftnbody) then
              mean_n(1:dimq) = mean_n(1:dimq) + qg(1:dimq)
            else
              do iiii = 1, dimq
                do jjjj = 1, dimq
                  covar_n(iiii, jjjj) = covar_n(iiii, jjjj) + (qg(iiii) - mean_n(iiii)) &
                                                            *(qg(jjjj) - mean_n(jjjj))
                end do
              end do
            end if
          else

            ! ---- derivatives of the Gram coordinates ----
            if (desc_forces_local) then
              dq(:, 1:dimq, 1:m) = 0.d0
              do p = 1, m
                dq(1:3, p, p) = d_ur_central(1:3, jp(p))
                ! partial cutoff products for the cutoff-derivative term
                cex(p) = 1.d0
                do q = 1, m
                  if (q /= p) cex(p) = cex(p)*fcp(q)
                end do
              end do
              do p = 1, m - 1
                do q = p + 1, m
                  ipos = gram_pair_pos(m, p, q)
                  cpq  = qg(ipos)
                  dq(1:3, ipos, p) = tmp_dxp(1:3, jp(q)+1)/(rp(p)*rp(q)) - tmp_dxp(1:3, jp(p)+1)*cpq/rp(p)**2
                  dq(1:3, ipos, q) = tmp_dxp(1:3, jp(p)+1)/(rp(p)*rp(q)) - tmp_dxp(1:3, jp(q)+1)*cpq/rp(q)**2
                end do
              end do
              ! R3: coordinate-derivative matrix for the dgemm contraction
              do t = 1, m
                do i = 1, dimq
                  DQmat(i, (t-1)*3+1:(t-1)*3+3) = dq(1:3, i, t)
                end do
              end do
            end if

            ! ---- cluster-level chemistry routing (permutation invariant) ----
            off = 0
            sgn = 1.d0
            if (ftnbody_chem_mode == 1 .or. ftnbody_chem_mode == 3) then
              select case (nord)
              case (2)
                call ftnbody_chem_route_2b(ftnbody_chem_mode, sa, sj(1), nsp, nrff, nchan, ic, off, sgn, chem_w)
              case (3)
                call ftnbody_chem_route_3b(ftnbody_chem_mode, sa, sj(1), sj(2), nsp, nrff, nchan, ic, off, sgn, chem_w)
              case (4)
                call ftnbody_chem_route_4b(ftnbody_chem_mode, sa, sj(1), sj(2), sj(3), nsp, nrff, nchan, ic, off, sgn, chem_w)
              case (5)
                call ftnbody_chem_route_5b(ftnbody_chem_mode, sa, sj(1), sj(2), sj(3), sj(4), nsp, nrff, nchan, ic, off, sgn, chem_w)
              end select
            end if

            ! ---- exact orbit average over the S_m neighbour permutations ----
            ! (batched: one dgemv for all permuted phases, one cos/sin, then the
            !  pi-sum is taken BEFORE the slot contraction and the channel
            !  accumulation, so the expensive per-(t,ix) work is done once per
            !  cluster instead of once per permutation)
            do ipm = 1, nper
              theta_all(1:nrff, ipm) = phase(1:nrff)
            end do
            call dgemv('T', dimq, nrff*nper, one, omegaP, dimq, qg, incx, one, theta_all, incy)
            cosA(1:nrff, 1:nper) = cos(theta_all(1:nrff, 1:nper))
            if (desc_forces_local) sinA(1:nrff, 1:nper) = sin(theta_all(1:nrff, 1:nper))

            ! joint geometry-chemistry symmetrization (mode 2): the same
            ! permutation acts on the geometry slots and on the species.
            ! Permutations with the same species ARRANGEMENT share their
            ! chemistry weights, so they are grouped and the orbit sums are
            ! shared (monospecies cluster -> a single group).  Modes 0/1/3 have
            ! permutation-independent channel weights -> always a single group.
            if (ftnbody_chem_mode == 2 .and. nsp > 1) then
              ngrp = 0
              do ipm = 1, nper
                k = 0
                do p = 1, m
                  k = k*(nsp + 1) + sj(gram_perm(p, ipm, nord))
                end do
                g = 0
                do i = 1, ngrp
                  if (grp_key(i) == k) then
                    g = i
                    exit
                  end if
                end do
                if (g == 0) then
                  ngrp = ngrp + 1
                  grp_key(ngrp) = k
                  g = ngrp
                  do irk = 1, nchan
                    wchem = ftnbody_chem_center(sa, irk, nord)
                    do p = 1, m
                      wchem = wchem*ftnbody_chem_neigh(sj(gram_perm(p, ipm, nord)), irk, p, nord)
                    end do
                    chem_wg(irk, g) = wchem
                  end do
                end if
                grp_of(ipm) = g
              end do
            else
              ngrp = 1
              grp_of(1:nper) = 1
              if (ftnbody_chem_mode == 2) then
                do irk = 1, nchan
                  wchem = ftnbody_chem_center(sa, irk, nord)
                  do p = 1, m
                    wchem = wchem*ftnbody_chem_neigh(sj(p), irk, p, nord)
                  end do
                  chem_wg(irk, 1) = wchem
                end do
              end if
            end if

            ! per-group orbit sums
            CEg(1:nrff, 1:ngrp) = 0.d0
            if (desc_forces_local) SJg(1:nrff, 1:dimq, 1:ngrp) = 0.d0
            do ipm = 1, nper
              g = grp_of(ipm)
              CEg(1:nrff, g) = CEg(1:nrff, g) + cosA(1:nrff, ipm)
              if (desc_forces_local) then
                do i = 1, dimq
                  SJg(1:nrff, i, g) = SJg(1:nrff, i, g) - sinA(1:nrff, ipm)*omegaPT(1:nrff, i, ipm)
                end do
              end if
            end do

            ! energy accumulation
            do g = 1, ngrp
              select case (ftnbody_chem_mode)
              case (0)
                wsc = winv*cprod
                tmp_real(1:nrff) = tmp_real(1:nrff) + wsc*CEg(1:nrff, g)
              case (2)
                do irk = 1, nchan
                  off = (irk - 1)*nrff
                  wsc = chem_wg(irk, g)*winv*cprod
                  tmp_real(off+1:off+nrff) = tmp_real(off+1:off+nrff) + wsc*CEg(1:nrff, g)
                end do
              case default   ! modes 1 (exact) and 3 (hash): one channel, weight sgn
                wsc = sgn*winv*cprod
                tmp_real(off+1:off+nrff) = tmp_real(off+1:off+nrff) + wsc*CEg(1:nrff, g)
              end select
            end do

            ! force accumulation (dgemm slot contraction once per cluster and group)
            if (desc_forces_local) then
              do g = 1, ngrp
                call dgemm('N', 'N', nrff, 3*m, dimq, one, SJg(1, 1, g), nrff, DQmat, MAXQ, zero, Wg, nrff)
                do t = 1, m
                  do ix = 1, 3
                    gF(1:nrff) = Wg(1:nrff, (t-1)*3+ix)*cprod &
                               + CEg(1:nrff, g)*cex(t)*dfcp(t)*tmp_dxp(ix, jp(t)+1)/rp(t)
                    select case (ftnbody_chem_mode)
                    case (0)
                      tmp_dreal(1:nrff, jp(t), ix) = tmp_dreal(1:nrff, jp(t), ix) + winv*gF(1:nrff)
                    case (2)
                      do irk = 1, nchan
                        off = (irk - 1)*nrff
                        tmp_dreal(off+1:off+nrff, jp(t), ix) = tmp_dreal(off+1:off+nrff, jp(t), ix) &
                                                             + chem_wg(irk, g)*winv*gF(1:nrff)
                      end do
                    case default
                      tmp_dreal(off+1:off+nrff, jp(t), ix) = tmp_dreal(off+1:off+nrff, jp(t), ix) &
                                                           + sgn*winv*gF(1:nrff)
                    end select
                  end do
                end do
              end do
            end if

          end if   ! init_mode / descriptor
        end if   ! okcut

        ! advance to the next m-combination of {1..max_neigh_local}
        k = m
        do
          if (k < 1) exit
          if (idx(k) < max_neigh_local - m + k) exit
          k = k - 1
        end do
        if (k < 1) exit odometer
        idx(k) = idx(k) + 1
        do p = k + 1, m
          idx(p) = idx(p-1) + 1
        end do
      end do odometer

      if (.not. init_mode_ftnbody) then
        tmp_real(:) = tmp_real(:)*norm
        if (desc_forces_local) then
          do k = 1, max_neigh_local
            do ix = 1, 3
              tmp_dreal(:, k, ix) = tmp_dreal(:, k, ix)*norm
            end do
          end do
        end if
      end if

      deallocate(fcut_all, dfcut_all, gF, Wg, chem_w)
      deallocate(omegaP, omegaPT, theta_all, cosA, sinA, CEg, SJg, chem_wg)
    end subroutine ftnbody_gram_generic

    ! ==================================================================
    ! ftnbody_model = 'cpip'  (docs/MLT5 sec. "Level B: cPIP --- compact
    ! mixed invariant representation")
    !
    ! Compact mixed invariant coordinates, versioned basis v1:
    !   n=2: (u1)                          n=3: (u1+u2, u1*u2, c12)
    !   n=4,5: (Su, Su2, Su3, Sc, Sc2, Sc3, P7, P8, P9) with
    !     P7 = sum_i u_i sum_{j/=i} c_ij
    !     P8 = sum_i u_i prod_{j<k, j,k/=i} c_jk
    !     P9 = sum_{i<j} u_i u_j c_ij
    ! Chemistry modes: 0 (blind: RFF input = P) and 2 (low rank: per rank r
    ! the colored moments C_r^u, C_r^u2, C_r^uc, C_r^c, C_r^uuc built from
    ! E^c_(sa,r) and the slot-independent neighbour embedding E_(s,r) are
    ! appended to P, and every rank owns one descriptor block).
    ! ------------------------------------------------------------------
    ! Dispatch wrapper: selects the per-order module arrays.
    ! ==================================================================
    subroutine ftnbody_cpip_order(nord, type_db_ja, max_neigh_local, i_type_db, r_central, tmp_dxp, &
                                  ur_central, d_ur_central, desc_forces_local)
      use module_ftnbody, only: ftnbody_omega_2b, ftnbody_omega_3b, ftnbody_omega_4b, ftnbody_omega_5b, &
                                ftnbody_phase_random_2b, ftnbody_phase_random_3b, &
                                ftnbody_phase_random_4b, ftnbody_phase_random_5b, &
                                r_cut_ft2b, r_cut_width_ft2b, r_cut_ft3b, r_cut_width_ft3b, &
                                r_cut_ft4b, r_cut_width_ft4b, r_cut_ft5b, r_cut_width_ft5b, &
                                count_2b, count_3b, count_4b, count_5b, &
                                mean_2b, mean_3b, mean_4b, mean_5b, &
                                covar_matrix_2b, covar_matrix_3b, covar_matrix_4b, covar_matrix_5b
      implicit none
      integer, intent(in) :: nord, type_db_ja, max_neigh_local
      integer, dimension(:), intent(in)  :: i_type_db
      real(kind_double), dimension(:), intent(in) :: r_central, ur_central
      real(kind_double), dimension(:,:), intent(in) :: tmp_dxp
      real(kind_double), dimension(:,:), intent(in) :: d_ur_central
      logical, intent(in) :: desc_forces_local

      select case (nord)
      case (2)
        call ftnbody_cpip_generic(2, type_db_ja, max_neigh_local, i_type_db, r_central, tmp_dxp, &
                                  ur_central, d_ur_central, desc_forces_local, &
                                  ftnbody_omega_2b, ftnbody_phase_random_2b, &
                                  r_cut_ft2b, r_cut_width_ft2b, count_2b, mean_2b, covar_matrix_2b)
      case (3)
        call ftnbody_cpip_generic(3, type_db_ja, max_neigh_local, i_type_db, r_central, tmp_dxp, &
                                  ur_central, d_ur_central, desc_forces_local, &
                                  ftnbody_omega_3b, ftnbody_phase_random_3b, &
                                  r_cut_ft3b, r_cut_width_ft3b, count_3b, mean_3b, covar_matrix_3b)
      case (4)
        call ftnbody_cpip_generic(4, type_db_ja, max_neigh_local, i_type_db, r_central, tmp_dxp, &
                                  ur_central, d_ur_central, desc_forces_local, &
                                  ftnbody_omega_4b, ftnbody_phase_random_4b, &
                                  r_cut_ft4b, r_cut_width_ft4b, count_4b, mean_4b, covar_matrix_4b)
      case (5)
        call ftnbody_cpip_generic(5, type_db_ja, max_neigh_local, i_type_db, r_central, tmp_dxp, &
                                  ur_central, d_ur_central, desc_forces_local, &
                                  ftnbody_omega_5b, ftnbody_phase_random_5b, &
                                  r_cut_ft5b, r_cut_width_ft5b, count_5b, mean_5b, covar_matrix_5b)
      end select
    end subroutine ftnbody_cpip_order

    ! ------------------------------------------------------------------
    ! Order-generic cPIP worker (energy descriptor block, analytic force
    ! block, and mean/covariance accumulation for the Fourier-length
    ! estimation passes).  The coordinates are already permutation
    ! invariant, so - unlike gramm - no orbit averaging is needed.
    ! ------------------------------------------------------------------
    subroutine ftnbody_cpip_generic(nord, type_db_ja, max_neigh_local, i_type_db, r_central, tmp_dxp, &
                                    ur_central, d_ur_central, desc_forces_local, &
                                    omega, phase, rcut_n, rcutw_n, cnt_n, mean_n, covar_n)
      use module_ftnbody, only: tmp_real, tmp_dreal, delta_rff, dim_rff, &
                                init_mode_ftnbody, ltmp_mean_ftnbody, &
                                ftnbody_chem_mode, ftnbody_n_channels, &
                                ftnbody_chem_center, ftnbody_chem_neigh, &
                                gram_pair_pos, CPIP_MG_N, CPIP_MC_N
      use module_neigh_local, only: r_cut_in, r_cut_width_in, r_cut_pair_in, fcut_rij_inout
      implicit none
      integer, intent(in) :: nord, type_db_ja, max_neigh_local
      integer, dimension(:), intent(in)  :: i_type_db
      real(kind_double), dimension(:), intent(in) :: r_central, ur_central
      real(kind_double), dimension(:,:), intent(in) :: tmp_dxp
      real(kind_double), dimension(:,:), intent(in) :: d_ur_central
      logical, intent(in) :: desc_forces_local
      real(kind_double), dimension(:,:), intent(in) :: omega   ! (dimq, nrff)
      real(kind_double), dimension(:),   intent(in) :: phase   ! (nrff)
      real(kind_double), intent(in)    :: rcut_n, rcutw_n
      real(kind_double), intent(inout) :: cnt_n
      real(kind_double), dimension(:),   intent(inout) :: mean_n
      real(kind_double), dimension(:,:), intent(inout) :: covar_n

      integer, parameter :: MAXM = 4, MAXX = 10, MAXH = 14
      real(kind_double), parameter :: one = 1.d0, zero = 0.d0
      integer :: m, dimx, mg, mc, dimq, nrff, nchan, sa
      integer :: idx(MAXM), jp(MAXM), sj(MAXM)
      real(kind_double) :: rp(MAXM), up(MAXM), fcp(MAXM), dfcp(MAXM), cex(MAXM)
      real(kind_double) :: Sang(MAXM), En(MAXM)
      real(kind_double) :: xg(MAXX)
      real(kind_double) :: dxprim(3, MAXX, MAXM)
      real(kind_double) :: Pg(MAXH), Cc(MAXH), hv(MAXH)
      ! transposed derivative blocks (dimx, moment): contiguous per moment
      real(kind_double) :: dPt(MAXX, MAXH), dCt(MAXX, MAXH)
      real(kind_double) :: DQmat(MAXX, 3*MAXM)
      real(kind_double) :: BgP(MAXH, 3*MAXM), BcC(MAXH, 3*MAXM)
      real(kind_double), dimension(:,:), allocatable :: omegaG, omegaC   ! (mg,nrff), (mc,nrff)
      real(kind_double), dimension(:),   allocatable :: theta0, theta, cosv, sinv, gF
      real(kind_double), dimension(:,:), allocatable :: Wg, Wc           ! (nrff, 3m)
      real(kind_double), dimension(:),   allocatable :: fcut_all, dfcut_all
      real(kind_double) :: cprod, cpq, Ec, local_r_cut_in, fcut, dfcut, norm
      integer :: p, q, t, i, ix, k, ipos, irk, off, ii, icol
      integer :: type_fcut_in, type_fcut_out, incx, incy
      logical :: okcut

      m    = nord - 1
      dimx = m + m*(m - 1)/2
      mg   = CPIP_MG_N(nord)
      mc   = 0
      if (ftnbody_chem_mode == 2) mc = CPIP_MC_N(nord)
      dimq = mg + mc
      nrff = dim_rff(nord)
      nchan = ftnbody_n_channels(nord)
      sa   = type_db_ja
      norm = 1.d0/sqrt(dble(nrff))*delta_rff(nord)
      incx = 1
      incy = 1

      tmp_real(:) = 0.d0
      if (desc_forces_local) tmp_dreal(:, :, :) = 0.d0
      if (max_neigh_local < m) return

      allocate(fcut_all(max_neigh_local), dfcut_all(max_neigh_local))
      ! contiguous copies of the geometric / colored omega row blocks
      allocate(omegaG(mg, nrff))
      omegaG(1:mg, 1:nrff) = omega(1:mg, 1:nrff)
      if (mc > 0) then
        allocate(omegaC(mc, nrff))
        omegaC(1:mc, 1:nrff) = omega(mg+1:dimq, 1:nrff)
      end if
      allocate(theta0(nrff), theta(nrff), cosv(nrff), sinv(nrff), gF(nrff))
      allocate(Wg(nrff, 3*MAXM), Wc(nrff, 3*MAXM))

      !TODOftnbody (same convention as the poly model)
      type_fcut_in  = 3
      type_fcut_out = 2
      do ii = 1, max_neigh_local
        ! Pair-specific inner cutoff: i_type_db(ii+1) due to assumed-shape remapping (index 1 = central)
        local_r_cut_in = r_cut_pair_in(type_db_ja, i_type_db(ii + 1))
        call fcut_rij_inout(type_fcut_in, type_fcut_out, r_central(ii), rcut_n, rcutw_n, &
                            local_r_cut_in, r_cut_width_in, desc_forces_local, fcut, dfcut)
        fcut_all(ii)  = fcut
        dfcut_all(ii) = dfcut
      end do

      ! enumerate every unordered cluster j_1 < j_2 < ... < j_m exactly once
      do p = 1, m
        idx(p) = p
      end do

      odometer: do
        okcut = .true.
        do p = 1, m
          if (fcut_all(idx(p)) == 0.d0) then
            okcut = .false.
            exit
          end if
        end do

        if (okcut) then
          ! ---- primitive coordinates of the cluster ----
          cprod = 1.d0
          do p = 1, m
            jp(p)  = idx(p)
            sj(p)  = i_type_db(jp(p) + 1)
            rp(p)  = r_central(jp(p))
            up(p)  = ur_central(jp(p))
            fcp(p) = fcut_all(jp(p))
            dfcp(p) = dfcut_all(jp(p))
            cprod  = cprod*fcp(p)
            xg(p)  = up(p)
          end do
          do p = 1, m - 1
            do q = p + 1, m
              cpq = dot_product(tmp_dxp(1:3, jp(p)+1), tmp_dxp(1:3, jp(q)+1))/(rp(p)*rp(q))
              xg(gram_pair_pos(m, p, q)) = cpq
            end do
          end do
          ! angular row sums S_p = sum_{q/=p} c_pq (used by P7 and C^uc)
          Sang(1:m) = 0.d0
          do p = 1, m - 1
            do q = p + 1, m
              cpq = xg(gram_pair_pos(m, p, q))
              Sang(p) = Sang(p) + cpq
              Sang(q) = Sang(q) + cpq
            end do
          end do

          call cpip_geom_moments(desc_forces_local)

          if (init_mode_ftnbody) then
            ! Fourier-length estimation passes: one sample per RFF input vector
            if (ftnbody_chem_mode == 2) then
              do irk = 1, nchan
                Ec = ftnbody_chem_center(sa, irk, nord)
                do p = 1, m
                  En(p) = ftnbody_chem_neigh(sj(p), irk, 1, nord)
                end do
                call cpip_colored_moments(.false.)
                hv(1:mg) = Pg(1:mg)
                hv(mg+1:dimq) = Cc(1:mc)
                call cpip_stats(hv)
              end do
            else
              hv(1:mg) = Pg(1:mg)
              call cpip_stats(hv)
            end if
          else

            ! ---- derivatives of the primitive coordinates ----
            if (desc_forces_local) then
              dxprim(:, 1:dimx, 1:m) = 0.d0
              do p = 1, m
                dxprim(1:3, p, p) = d_ur_central(1:3, jp(p))
                ! partial cutoff products for the cutoff-derivative term
                cex(p) = 1.d0
                do q = 1, m
                  if (q /= p) cex(p) = cex(p)*fcp(q)
                end do
              end do
              do p = 1, m - 1
                do q = p + 1, m
                  ipos = gram_pair_pos(m, p, q)
                  cpq  = xg(ipos)
                  dxprim(1:3, ipos, p) = tmp_dxp(1:3, jp(q)+1)/(rp(p)*rp(q)) - tmp_dxp(1:3, jp(p)+1)*cpq/rp(p)**2
                  dxprim(1:3, ipos, q) = tmp_dxp(1:3, jp(p)+1)/(rp(p)*rp(q)) - tmp_dxp(1:3, jp(q)+1)*cpq/rp(q)**2
                end do
              end do
              do t = 1, m
                do i = 1, dimx
                  DQmat(i, (t-1)*3+1:(t-1)*3+3) = dxprim(1:3, i, t)
                end do
              end do
              ! geometric part of the slot contraction (rank independent):
              ! Wg = omegaG^T * (dP * DQ)
              call dgemm('T', 'N', mg, 3*m, dimx, one, dPt, MAXX, DQmat, MAXX, zero, BgP, MAXH)
              call dgemm('T', 'N', nrff, 3*m, mg, one, omegaG, mg, BgP, MAXH, zero, Wg, nrff)
            end if

            ! geometric part of the Fourier phase (rank independent)
            theta0(1:nrff) = phase(1:nrff)
            call dgemv('T', mg, nrff, one, omegaG, mg, Pg, incx, one, theta0, incy)

            if (ftnbody_chem_mode == 0) then
              cosv(:) = cos(theta0(:))
              tmp_real(1:nrff) = tmp_real(1:nrff) + cosv(1:nrff)*cprod
              if (desc_forces_local) then
                sinv(:) = sin(theta0(:))
                do t = 1, m
                  do ix = 1, 3
                    icol = (t-1)*3 + ix
                    gF(1:nrff) = -sinv(1:nrff)*Wg(1:nrff, icol)*cprod &
                               + cosv(1:nrff)*cex(t)*dfcp(t)*tmp_dxp(ix, jp(t)+1)/rp(t)
                    tmp_dreal(1:nrff, jp(t), ix) = tmp_dreal(1:nrff, jp(t), ix) + gF(1:nrff)
                  end do
                end do
              end if
            else
              ! ---- chem mode 2: one descriptor block per rank, colored moments
              !      appended to the RFF input (joint chemistry-geometry coupling)
              do irk = 1, nchan
                Ec = ftnbody_chem_center(sa, irk, nord)
                do p = 1, m
                  En(p) = ftnbody_chem_neigh(sj(p), irk, 1, nord)
                end do
                call cpip_colored_moments(desc_forces_local)
                theta(1:nrff) = theta0(1:nrff)
                call dgemv('T', mc, nrff, one, omegaC, mc, Cc, incx, one, theta, incy)
                cosv(:) = cos(theta(:))
                off = (irk - 1)*nrff
                tmp_real(off+1:off+nrff) = tmp_real(off+1:off+nrff) + cosv(1:nrff)*cprod
                if (desc_forces_local) then
                  sinv(:) = sin(theta(:))
                  call dgemm('T', 'N', mc, 3*m, dimx, one, dCt, MAXX, DQmat, MAXX, zero, BcC, MAXH)
                  call dgemm('T', 'N', nrff, 3*m, mc, one, omegaC, mc, BcC, MAXH, zero, Wc, nrff)
                  do t = 1, m
                    do ix = 1, 3
                      icol = (t-1)*3 + ix
                      gF(1:nrff) = -sinv(1:nrff)*(Wg(1:nrff, icol) + Wc(1:nrff, icol))*cprod &
                                 + cosv(1:nrff)*cex(t)*dfcp(t)*tmp_dxp(ix, jp(t)+1)/rp(t)
                      tmp_dreal(off+1:off+nrff, jp(t), ix) = tmp_dreal(off+1:off+nrff, jp(t), ix) + gF(1:nrff)
                    end do
                  end do
                end if
              end do
            end if

          end if   ! init_mode / descriptor
        end if   ! okcut

        ! advance to the next m-combination of {1..max_neigh_local}
        k = m
        do
          if (k < 1) exit
          if (idx(k) < max_neigh_local - m + k) exit
          k = k - 1
        end do
        if (k < 1) exit odometer
        idx(k) = idx(k) + 1
        do p = k + 1, m
          idx(p) = idx(p-1) + 1
        end do
      end do odometer

      if (.not. init_mode_ftnbody) then
        tmp_real(:) = tmp_real(:)*norm
        if (desc_forces_local) then
          do k = 1, max_neigh_local
            do ix = 1, 3
              tmp_dreal(:, k, ix) = tmp_dreal(:, k, ix)*norm
            end do
          end do
        end if
      end if

      deallocate(fcut_all, dfcut_all, omegaG, theta0, theta, cosv, sinv, gF, Wg, Wc)
      if (allocated(omegaC)) deallocate(omegaC)

    contains

      ! geometric moments P (versioned compact basis v1) and, on request,
      ! their derivatives w.r.t. the primitive coordinates (dPt(dimx, moment))
      subroutine cpip_geom_moments(want_deriv)
        logical, intent(in) :: want_deriv
        integer :: lp, lq, la, lb, lia, lib, lpos
        real(kind_double) :: cval, pr, ds
        Pg(1:mg) = 0.d0
        if (want_deriv) dPt(1:dimx, 1:mg) = 0.d0
        if (m == 1) then
          Pg(1) = up(1)
          if (want_deriv) dPt(1, 1) = 1.d0
          return
        end if
        if (m == 2) then
          ! complete 3-body invariant set (identical content to the PIP level)
          Pg(1) = up(1) + up(2)
          Pg(2) = up(1)*up(2)
          Pg(3) = xg(3)
          if (want_deriv) then
            dPt(1, 1) = 1.d0
            dPt(2, 1) = 1.d0
            dPt(1, 2) = up(2)
            dPt(2, 2) = up(1)
            dPt(3, 3) = 1.d0
          end if
          return
        end if
        ! m >= 3: radial / angular power sums
        do lp = 1, m
          Pg(1) = Pg(1) + up(lp)
          Pg(2) = Pg(2) + up(lp)**2
          Pg(3) = Pg(3) + up(lp)**3
          if (want_deriv) then
            dPt(lp, 1) = 1.d0
            dPt(lp, 2) = 2.d0*up(lp)
            dPt(lp, 3) = 3.d0*up(lp)**2
          end if
        end do
        do lpos = m + 1, dimx
          cval = xg(lpos)
          Pg(4) = Pg(4) + cval
          Pg(5) = Pg(5) + cval**2
          Pg(6) = Pg(6) + cval**3
          if (want_deriv) then
            dPt(lpos, 4) = 1.d0
            dPt(lpos, 5) = 2.d0*cval
            dPt(lpos, 6) = 3.d0*cval**2
          end if
        end do
        ! P7 = sum_p u_p S_p ; P9 = sum_{p<q} u_p u_q c_pq
        do lp = 1, m
          Pg(7) = Pg(7) + up(lp)*Sang(lp)
          if (want_deriv) dPt(lp, 7) = Sang(lp)
        end do
        do lp = 1, m - 1
          do lq = lp + 1, m
            lpos = gram_pair_pos(m, lp, lq)
            cval = xg(lpos)
            Pg(9) = Pg(9) + up(lp)*up(lq)*cval
            if (want_deriv) then
              dPt(lpos, 7) = up(lp) + up(lq)
              dPt(lpos, 9) = up(lp)*up(lq)
              dPt(lp, 9) = dPt(lp, 9) + up(lq)*cval
              dPt(lq, 9) = dPt(lq, 9) + up(lp)*cval
            end if
          end do
        end do
        ! P8 = sum_p u_p prod_{a<b, a,b/=p} c_ab
        do lp = 1, m
          pr = 1.d0
          do la = 1, m - 1
            if (la == lp) cycle
            do lb = la + 1, m
              if (lb == lp) cycle
              pr = pr*xg(gram_pair_pos(m, la, lb))
            end do
          end do
          Pg(8) = Pg(8) + up(lp)*pr
          if (want_deriv) dPt(lp, 8) = pr
        end do
        if (want_deriv) then
          do la = 1, m - 1
            do lb = la + 1, m
              lpos = gram_pair_pos(m, la, lb)
              ds = 0.d0
              do lp = 1, m
                if (lp == la .or. lp == lb) cycle
                pr = 1.d0
                do lia = 1, m - 1
                  if (lia == lp) cycle
                  do lib = lia + 1, m
                    if (lib == lp) cycle
                    if (lia == la .and. lib == lb) cycle
                    pr = pr*xg(gram_pair_pos(m, lia, lib))
                  end do
                end do
                ds = ds + up(lp)*pr
              end do
              dPt(lpos, 8) = ds
            end do
          end do
        end if
      end subroutine cpip_geom_moments

      ! colored chemistry-geometry moments of rank irk (Ec, En already set):
      !   C^u = Ec sum_p En_p u_p            C^u2 = Ec sum_p En_p u_p^2
      !   C^uc = Ec sum_p En_p u_p S_p       C^c = Ec sum_{p<q} En_p En_q c_pq
      !   C^uuc = Ec sum_{p<q} En_p En_q u_p u_q c_pq
      ! (m=1 keeps only C^u, C^u2); derivatives in dCt(dimx, moment) on request
      subroutine cpip_colored_moments(want_deriv)
        logical, intent(in) :: want_deriv
        integer :: lp, lq, lpos
        real(kind_double) :: cval, ee, epq
        Cc(1:mc) = 0.d0
        if (want_deriv) dCt(1:dimx, 1:mc) = 0.d0
        if (m == 1) then
          ee = Ec*En(1)
          Cc(1) = ee*up(1)
          Cc(2) = ee*up(1)**2
          if (want_deriv) then
            dCt(1, 1) = ee
            dCt(1, 2) = 2.d0*ee*up(1)
          end if
          return
        end if
        do lp = 1, m
          ee = Ec*En(lp)
          Cc(1) = Cc(1) + ee*up(lp)
          Cc(2) = Cc(2) + ee*up(lp)**2
          Cc(3) = Cc(3) + ee*up(lp)*Sang(lp)
          if (want_deriv) then
            dCt(lp, 1) = ee
            dCt(lp, 2) = 2.d0*ee*up(lp)
            dCt(lp, 3) = ee*Sang(lp)
          end if
        end do
        do lp = 1, m - 1
          do lq = lp + 1, m
            lpos = gram_pair_pos(m, lp, lq)
            cval = xg(lpos)
            epq = Ec*En(lp)*En(lq)
            Cc(4) = Cc(4) + epq*cval
            Cc(5) = Cc(5) + epq*up(lp)*up(lq)*cval
            if (want_deriv) then
              dCt(lpos, 3) = Ec*(En(lp)*up(lp) + En(lq)*up(lq))
              dCt(lpos, 4) = epq
              dCt(lpos, 5) = epq*up(lp)*up(lq)
              dCt(lp, 5) = dCt(lp, 5) + epq*up(lq)*cval
              dCt(lq, 5) = dCt(lq, 5) + epq*up(lp)*cval
            end if
          end do
        end do
      end subroutine cpip_colored_moments

      ! mean / covariance accumulation of one RFF input sample
      subroutine cpip_stats(hs)
        real(kind_double), intent(in) :: hs(:)
        integer :: li, lj
        cnt_n = cnt_n + 1
        if (ltmp_mean_ftnbody) then
          mean_n(1:dimq) = mean_n(1:dimq) + hs(1:dimq)
        else
          do li = 1, dimq
            do lj = 1, dimq
              covar_n(li, lj) = covar_n(li, lj) + (hs(li) - mean_n(li))*(hs(lj) - mean_n(lj))
            end do
          end do
        end if
      end subroutine cpip_stats

    end subroutine ftnbody_cpip_generic

    ! ==================================================================
    ! ftnbody_model = 'spip'  (docs/MLT5 sec. "Level C: sPIP --- systematic
    ! orbit-polynomial representation")
    !
    ! Systematic S_{n-1} monomial-orbit invariant coordinates up to
    ! spip_degree_nbody, evaluated from the tables built once in
    ! init_ftnbody_spip_orbits.  Chemistry modes:
    !   0: RFF input = orbit vector P (collapsed member list, no group op)
    !   1: exact colored orbit PIPs Q_{alpha,gamma} (restricted orbit sum over
    !      the permutations matching the sorted species representative), one
    !      descriptor block per exact channel (routing via chem_route_*b)
    !   2: low-rank colored orbit PIPs Q_{alpha,r} (joint symmetrization:
    !      the same permutation acts on the geometry and on the role
    !      embeddings), one block per rank; permutations sharing a species
    !      arrangement are grouped, monospecies clusters use the collapsed
    !      member list directly
    !   3: hashed exact channels (same coordinates as mode 1, hashed routing)
    ! ------------------------------------------------------------------
    ! Dispatch wrapper: selects the per-order module arrays.
    ! ==================================================================
    subroutine ftnbody_spip_order(nord, type_db_ja, max_neigh_local, i_type_db, r_central, tmp_dxp, &
                                  ur_central, d_ur_central, desc_forces_local)
      use module_ftnbody, only: ftnbody_omega_2b, ftnbody_omega_3b, ftnbody_omega_4b, ftnbody_omega_5b, &
                                ftnbody_phase_random_2b, ftnbody_phase_random_3b, &
                                ftnbody_phase_random_4b, ftnbody_phase_random_5b, &
                                r_cut_ft2b, r_cut_width_ft2b, r_cut_ft3b, r_cut_width_ft3b, &
                                r_cut_ft4b, r_cut_width_ft4b, r_cut_ft5b, r_cut_width_ft5b, &
                                count_2b, count_3b, count_4b, count_5b, &
                                mean_2b, mean_3b, mean_4b, mean_5b, &
                                covar_matrix_2b, covar_matrix_3b, covar_matrix_4b, covar_matrix_5b
      implicit none
      integer, intent(in) :: nord, type_db_ja, max_neigh_local
      integer, dimension(:), intent(in)  :: i_type_db
      real(kind_double), dimension(:), intent(in) :: r_central, ur_central
      real(kind_double), dimension(:,:), intent(in) :: tmp_dxp
      real(kind_double), dimension(:,:), intent(in) :: d_ur_central
      logical, intent(in) :: desc_forces_local

      select case (nord)
      case (2)
        call ftnbody_spip_generic(2, type_db_ja, max_neigh_local, i_type_db, r_central, tmp_dxp, &
                                  ur_central, d_ur_central, desc_forces_local, &
                                  ftnbody_omega_2b, ftnbody_phase_random_2b, &
                                  r_cut_ft2b, r_cut_width_ft2b, count_2b, mean_2b, covar_matrix_2b)
      case (3)
        call ftnbody_spip_generic(3, type_db_ja, max_neigh_local, i_type_db, r_central, tmp_dxp, &
                                  ur_central, d_ur_central, desc_forces_local, &
                                  ftnbody_omega_3b, ftnbody_phase_random_3b, &
                                  r_cut_ft3b, r_cut_width_ft3b, count_3b, mean_3b, covar_matrix_3b)
      case (4)
        call ftnbody_spip_generic(4, type_db_ja, max_neigh_local, i_type_db, r_central, tmp_dxp, &
                                  ur_central, d_ur_central, desc_forces_local, &
                                  ftnbody_omega_4b, ftnbody_phase_random_4b, &
                                  r_cut_ft4b, r_cut_width_ft4b, count_4b, mean_4b, covar_matrix_4b)
      case (5)
        call ftnbody_spip_generic(5, type_db_ja, max_neigh_local, i_type_db, r_central, tmp_dxp, &
                                  ur_central, d_ur_central, desc_forces_local, &
                                  ftnbody_omega_5b, ftnbody_phase_random_5b, &
                                  r_cut_ft5b, r_cut_width_ft5b, count_5b, mean_5b, covar_matrix_5b)
      end select
    end subroutine ftnbody_spip_order

    ! ------------------------------------------------------------------
    ! Order-generic sPIP worker (energy descriptor block, analytic force
    ! block, and mean/covariance accumulation for the Fourier-length
    ! estimation passes).
    ! ------------------------------------------------------------------
    subroutine ftnbody_spip_generic(nord, type_db_ja, max_neigh_local, i_type_db, r_central, tmp_dxp, &
                                    ur_central, d_ur_central, desc_forces_local, &
                                    omega, phase, rcut_n, rcutw_n, cnt_n, mean_n, covar_n)
      use module_ftnbody, only: tmp_real, tmp_dreal, delta_rff, dim_rff, &
                                init_mode_ftnbody, ltmp_mean_ftnbody, &
                                ftnbody_chem_mode, ftnbody_n_channels, &
                                ftnbody_chem_center, ftnbody_chem_neigh, &
                                gram_nperm, gram_perm, gram_qmap, gram_pair_pos, &
                                spip_tab, spip_degree_n
      use module_chemical_species, only: fix_no_of_elements
      use module_neigh_local, only: r_cut_in, r_cut_width_in, r_cut_pair_in, fcut_rij_inout
      implicit none
      integer, intent(in) :: nord, type_db_ja, max_neigh_local
      integer, dimension(:), intent(in)  :: i_type_db
      real(kind_double), dimension(:), intent(in) :: r_central, ur_central
      real(kind_double), dimension(:,:), intent(in) :: tmp_dxp
      real(kind_double), dimension(:,:), intent(in) :: d_ur_central
      logical, intent(in) :: desc_forces_local
      real(kind_double), dimension(:,:), intent(in) :: omega   ! (norb, nrff)
      real(kind_double), dimension(:),   intent(in) :: phase   ! (nrff)
      real(kind_double), intent(in)    :: rcut_n, rcutw_n
      real(kind_double), intent(inout) :: cnt_n
      real(kind_double), dimension(:),   intent(inout) :: mean_n
      real(kind_double), dimension(:,:), intent(inout) :: covar_n

      integer, parameter :: MAXM = 4, MAXX = 10, MAXPER = 24, MAXDEG = 10
      real(kind_double), parameter :: one = 1.d0, zero = 0.d0
      integer :: m, dimx, ndeg, norb, nrff, nper, nchan, nsp, sa, ngrp
      integer :: idx(MAXM), jp(MAXM), sj(MAXM), sbar(MAXM)
      integer :: grp_of(MAXPER), grp_key(MAXPER)
      real(kind_double) :: rp(MAXM), up(MAXM), fcp(MAXM), dfcp(MAXM), cex(MAXM)
      real(kind_double) :: xg(MAXX)
      real(kind_double) :: dxprim(3, MAXX, MAXM)
      real(kind_double) :: DQmat(MAXX, 3*MAXM)
      real(kind_double) :: powx(0:MAXDEG, MAXX)
      real(kind_double) :: chem_dummy(1)
      integer :: eimg(MAXX)
      real(kind_double), dimension(:),     allocatable :: Qv, theta, cosv, sinv, gF
      real(kind_double), dimension(:,:),   allocatable :: dQt          ! (dimx, norb)
      real(kind_double), dimension(:,:),   allocatable :: Gv           ! (norb, nper)
      real(kind_double), dimension(:,:,:), allocatable :: dGt          ! (dimx, norb, nper)
      real(kind_double), dimension(:,:),   allocatable :: chem_wg      ! (nchan, nper)
      real(kind_double), dimension(:,:),   allocatable :: Bmat         ! (norb, 3m)
      real(kind_double), dimension(:,:),   allocatable :: Wmat         ! (nrff, 3m)
      real(kind_double), dimension(:,:,:), allocatable :: Wgm          ! (nrff, 3m, nper)
      real(kind_double), dimension(:),     allocatable :: fcut_all, dfcut_all
      real(kind_double) :: cprod, cpq, wchem, sgn, local_r_cut_in, fcut, dfcut, norm, winv
      integer :: p, q, t, i, ix, k, ipos, irk, off, ii, ic, ipm, g, iorb, imono, icol
      integer :: type_fcut_in, type_fcut_out, incx, incy
      logical :: okcut, match

      m    = nord - 1
      dimx = m + m*(m - 1)/2
      ndeg = spip_degree_n(nord)
      norb = spip_tab(nord)%norb
      nrff = dim_rff(nord)
      nper = gram_nperm(nord)
      winv = 1.d0/dble(nper)
      nsp  = fix_no_of_elements
      nchan = ftnbody_n_channels(nord)
      sa   = type_db_ja
      norm = 1.d0/sqrt(dble(nrff))*delta_rff(nord)
      incx = 1
      incy = 1

      tmp_real(:) = 0.d0
      if (desc_forces_local) tmp_dreal(:, :, :) = 0.d0
      if (max_neigh_local < m) return

      allocate(fcut_all(max_neigh_local), dfcut_all(max_neigh_local))
      allocate(Qv(norb), theta(nrff), cosv(nrff), sinv(nrff), gF(nrff))
      allocate(dQt(dimx, norb), Bmat(norb, 3*MAXM), Wmat(nrff, 3*MAXM))
      if (ftnbody_chem_mode == 2) then
        allocate(Gv(norb, nper), chem_wg(max(nchan, 1), nper))
        if (desc_forces_local) allocate(dGt(dimx, norb, nper), Wgm(nrff, 3*MAXM, nper))
      end if

      !TODOftnbody (same convention as the poly model)
      type_fcut_in  = 3
      type_fcut_out = 2
      do ii = 1, max_neigh_local
        ! Pair-specific inner cutoff: i_type_db(ii+1) due to assumed-shape remapping (index 1 = central)
        local_r_cut_in = r_cut_pair_in(type_db_ja, i_type_db(ii + 1))
        call fcut_rij_inout(type_fcut_in, type_fcut_out, r_central(ii), rcut_n, rcutw_n, &
                            local_r_cut_in, r_cut_width_in, desc_forces_local, fcut, dfcut)
        fcut_all(ii)  = fcut
        dfcut_all(ii) = dfcut
      end do

      ! enumerate every unordered cluster j_1 < j_2 < ... < j_m exactly once
      do p = 1, m
        idx(p) = p
      end do

      odometer: do
        okcut = .true.
        do p = 1, m
          if (fcut_all(idx(p)) == 0.d0) then
            okcut = .false.
            exit
          end if
        end do

        if (okcut) then
          ! ---- primitive coordinates of the cluster ----
          cprod = 1.d0
          do p = 1, m
            jp(p)  = idx(p)
            sj(p)  = i_type_db(jp(p) + 1)
            rp(p)  = r_central(jp(p))
            up(p)  = ur_central(jp(p))
            fcp(p) = fcut_all(jp(p))
            dfcp(p) = dfcut_all(jp(p))
            cprod  = cprod*fcp(p)
            xg(p)  = up(p)
          end do
          do p = 1, m - 1
            do q = p + 1, m
              cpq = dot_product(tmp_dxp(1:3, jp(p)+1), tmp_dxp(1:3, jp(q)+1))/(rp(p)*rp(q))
              xg(gram_pair_pos(m, p, q)) = cpq
            end do
          end do
          ! power table x_i^d, d = 0..ndeg (monomials are evaluated from it)
          do i = 1, dimx
            powx(0, i) = 1.d0
            do k = 1, ndeg
              powx(k, i) = powx(k-1, i)*xg(i)
            end do
          end do

          ! ---- derivatives of the primitive coordinates ----
          if (desc_forces_local) then
            dxprim(:, 1:dimx, 1:m) = 0.d0
            do p = 1, m
              dxprim(1:3, p, p) = d_ur_central(1:3, jp(p))
              ! partial cutoff products for the cutoff-derivative term
              cex(p) = 1.d0
              do q = 1, m
                if (q /= p) cex(p) = cex(p)*fcp(q)
              end do
            end do
            do p = 1, m - 1
              do q = p + 1, m
                ipos = gram_pair_pos(m, p, q)
                cpq  = xg(ipos)
                dxprim(1:3, ipos, p) = tmp_dxp(1:3, jp(q)+1)/(rp(p)*rp(q)) - tmp_dxp(1:3, jp(p)+1)*cpq/rp(p)**2
                dxprim(1:3, ipos, q) = tmp_dxp(1:3, jp(p)+1)/(rp(p)*rp(q)) - tmp_dxp(1:3, jp(q)+1)*cpq/rp(q)**2
              end do
            end do
            do t = 1, m
              do i = 1, dimx
                DQmat(i, (t-1)*3+1:(t-1)*3+3) = dxprim(1:3, i, t)
              end do
            end do
          end if

          select case (ftnbody_chem_mode)
          case (0)
            ! ---- chemistry blind: full orbit sums from the collapsed lists ----
            Qv(1:norb) = 0.d0
            if (desc_forces_local) dQt(1:dimx, 1:norb) = 0.d0
            do iorb = 1, norb
              do imono = spip_tab(nord)%orb_ptr(iorb), spip_tab(nord)%orb_ptr(iorb+1) - 1
                call spip_accum_mono(spip_tab(nord)%mono_expo(1:dimx, imono), &
                                     spip_tab(nord)%mono_coef(imono), desc_forces_local, &
                                     Qv(iorb), dQt(1:dimx, iorb))
              end do
            end do
            if (init_mode_ftnbody) then
              call spip_stats(Qv)
            else
              call spip_channel_accum(0, 1.d0)
            end if

          case (1, 3)
            ! ---- exact / hashed colored orbit PIPs: restricted orbit sum over
            !      the permutations matching the sorted species representative
            sbar(1:m) = sj(1:m)
            call isort_m(sbar, m)
            sgn = 1.d0
            ic = 1
            off = 0
            select case (nord)
            case (2)
              call ftnbody_chem_route_2b(ftnbody_chem_mode, sa, sj(1), nsp, nrff, nchan, ic, off, sgn, chem_dummy)
            case (3)
              call ftnbody_chem_route_3b(ftnbody_chem_mode, sa, sj(1), sj(2), nsp, nrff, nchan, ic, off, sgn, chem_dummy)
            case (4)
              call ftnbody_chem_route_4b(ftnbody_chem_mode, sa, sj(1), sj(2), sj(3), nsp, nrff, nchan, ic, off, sgn, chem_dummy)
            case (5)
              call ftnbody_chem_route_5b(ftnbody_chem_mode, sa, sj(1), sj(2), sj(3), sj(4), nsp, nrff, nchan, ic, off, sgn, chem_dummy)
            end select
            Qv(1:norb) = 0.d0
            if (desc_forces_local) dQt(1:dimx, 1:norb) = 0.d0
            do ipm = 1, nper
              match = .true.
              do p = 1, m
                if (sj(gram_perm(p, ipm, nord)) /= sbar(p)) then
                  match = .false.
                  exit
                end if
              end do
              if (.not. match) cycle
              do iorb = 1, norb
                do i = 1, dimx
                  eimg(gram_qmap(i, ipm, nord)) = spip_tab(nord)%rep_expo(i, iorb)
                end do
                call spip_accum_mono(eimg(1:dimx), winv, desc_forces_local, &
                                     Qv(iorb), dQt(1:dimx, iorb))
              end do
            end do
            if (init_mode_ftnbody) then
              call spip_stats(Qv)
            else
              call spip_channel_accum(off, sgn)
            end if

          case (2)
            ! ---- low-rank colored orbit PIPs: joint geometry-chemistry
            !      symmetrization; permutations sharing a species arrangement
            !      are grouped (monospecies cluster -> single group)
            if (nsp > 1) then
              ngrp = 0
              do ipm = 1, nper
                k = 0
                do p = 1, m
                  k = k*(nsp + 1) + sj(gram_perm(p, ipm, nord))
                end do
                g = 0
                do i = 1, ngrp
                  if (grp_key(i) == k) then
                    g = i
                    exit
                  end if
                end do
                if (g == 0) then
                  ngrp = ngrp + 1
                  grp_key(ngrp) = k
                  g = ngrp
                  do irk = 1, nchan
                    wchem = ftnbody_chem_center(sa, irk, nord)
                    do p = 1, m
                      wchem = wchem*ftnbody_chem_neigh(sj(gram_perm(p, ipm, nord)), irk, p, nord)
                    end do
                    chem_wg(irk, g) = wchem
                  end do
                end if
                grp_of(ipm) = g
              end do
            else
              ngrp = 1
              grp_of(1:nper) = 1
              do irk = 1, nchan
                wchem = ftnbody_chem_center(sa, irk, nord)
                do p = 1, m
                  wchem = wchem*ftnbody_chem_neigh(sj(p), irk, p, nord)
                end do
                chem_wg(irk, 1) = wchem
              end do
            end if

            ! restricted orbit sums per group: G(alpha,g) = (1/m!) sum_{pi in g} M_alpha(Pi_pi x)
            Gv(1:norb, 1:ngrp) = 0.d0
            if (desc_forces_local) dGt(1:dimx, 1:norb, 1:ngrp) = 0.d0
            if (ngrp == 1) then
              ! single species arrangement: the group sum is the full Reynolds
              ! average -> use the collapsed member list (no permutation loop)
              do iorb = 1, norb
                do imono = spip_tab(nord)%orb_ptr(iorb), spip_tab(nord)%orb_ptr(iorb+1) - 1
                  call spip_accum_mono(spip_tab(nord)%mono_expo(1:dimx, imono), &
                                       spip_tab(nord)%mono_coef(imono), desc_forces_local, &
                                       Gv(iorb, 1), dGt(1:dimx, iorb, 1))
                end do
              end do
            else
              do ipm = 1, nper
                g = grp_of(ipm)
                do iorb = 1, norb
                  do i = 1, dimx
                    eimg(gram_qmap(i, ipm, nord)) = spip_tab(nord)%rep_expo(i, iorb)
                  end do
                  call spip_accum_mono(eimg(1:dimx), winv, desc_forces_local, &
                                       Gv(iorb, g), dGt(1:dimx, iorb, g))
                end do
              end do
            end if

            if (init_mode_ftnbody) then
              do irk = 1, nchan
                Qv(1:norb) = 0.d0
                do g = 1, ngrp
                  Qv(1:norb) = Qv(1:norb) + chem_wg(irk, g)*Gv(1:norb, g)
                end do
                call spip_stats(Qv)
              end do
            else
              ! group slot contractions, once per cluster
              if (desc_forces_local) then
                do g = 1, ngrp
                  call dgemm('T', 'N', norb, 3*m, dimx, one, dGt(1, 1, g), dimx, DQmat, MAXX, zero, Bmat, norb)
                  call dgemm('T', 'N', nrff, 3*m, norb, one, omega, norb, Bmat, norb, zero, Wgm(1, 1, g), nrff)
                end do
              end if
              do irk = 1, nchan
                Qv(1:norb) = 0.d0
                do g = 1, ngrp
                  Qv(1:norb) = Qv(1:norb) + chem_wg(irk, g)*Gv(1:norb, g)
                end do
                theta(1:nrff) = phase(1:nrff)
                call dgemv('T', norb, nrff, one, omega, norb, Qv, incx, one, theta, incy)
                cosv(:) = cos(theta(:))
                off = (irk - 1)*nrff
                tmp_real(off+1:off+nrff) = tmp_real(off+1:off+nrff) + cosv(1:nrff)*cprod
                if (desc_forces_local) then
                  sinv(:) = sin(theta(:))
                  Wmat(1:nrff, 1:3*m) = 0.d0
                  do g = 1, ngrp
                    Wmat(1:nrff, 1:3*m) = Wmat(1:nrff, 1:3*m) + chem_wg(irk, g)*Wgm(1:nrff, 1:3*m, g)
                  end do
                  do t = 1, m
                    do ix = 1, 3
                      icol = (t-1)*3 + ix
                      gF(1:nrff) = -sinv(1:nrff)*Wmat(1:nrff, icol)*cprod &
                                 + cosv(1:nrff)*cex(t)*dfcp(t)*tmp_dxp(ix, jp(t)+1)/rp(t)
                      tmp_dreal(off+1:off+nrff, jp(t), ix) = tmp_dreal(off+1:off+nrff, jp(t), ix) + gF(1:nrff)
                    end do
                  end do
                end if
              end do
            end if
          end select

        end if   ! okcut

        ! advance to the next m-combination of {1..max_neigh_local}
        k = m
        do
          if (k < 1) exit
          if (idx(k) < max_neigh_local - m + k) exit
          k = k - 1
        end do
        if (k < 1) exit odometer
        idx(k) = idx(k) + 1
        do p = k + 1, m
          idx(p) = idx(p-1) + 1
        end do
      end do odometer

      if (.not. init_mode_ftnbody) then
        tmp_real(:) = tmp_real(:)*norm
        if (desc_forces_local) then
          do k = 1, max_neigh_local
            do ix = 1, 3
              tmp_dreal(:, k, ix) = tmp_dreal(:, k, ix)*norm
            end do
          end do
        end if
      end if

      deallocate(fcut_all, dfcut_all, Qv, theta, cosv, sinv, gF, dQt, Bmat, Wmat)
      if (allocated(Gv))      deallocate(Gv)
      if (allocated(dGt))     deallocate(dGt)
      if (allocated(chem_wg)) deallocate(chem_wg)
      if (allocated(Wgm))     deallocate(Wgm)

    contains

      ! accumulate one monomial coef * prod_i x_i^e_i (values from the host
      ! powx table) into qa, and its primitive-coordinate gradient into dqa
      subroutine spip_accum_mono(evec, coef, want_deriv, qa, dqa)
        integer, intent(in) :: evec(:)
        real(kind_double), intent(in) :: coef
        logical, intent(in) :: want_deriv
        real(kind_double), intent(inout) :: qa
        real(kind_double), intent(inout) :: dqa(:)
        integer :: nz(MAXX), nnz, li, lj, lv
        real(kind_double) :: v, dv
        nnz = 0
        v = coef
        do li = 1, dimx
          if (evec(li) > 0) then
            nnz = nnz + 1
            nz(nnz) = li
            v = v*powx(evec(li), li)
          end if
        end do
        qa = qa + v
        if (.not. want_deriv) return
        do li = 1, nnz
          lv = nz(li)
          dv = coef*dble(evec(lv))*powx(evec(lv)-1, lv)
          do lj = 1, nnz
            if (lj == li) cycle
            dv = dv*powx(evec(nz(lj)), nz(lj))
          end do
          dqa(lv) = dqa(lv) + dv
        end do
      end subroutine spip_accum_mono

      ! energy + force accumulation of the current Qv/dQt into the descriptor
      ! block starting at loff, with channel weight lsgn (modes 0, 1 and 3)
      subroutine spip_channel_accum(loff, lsgn)
        integer, intent(in) :: loff
        real(kind_double), intent(in) :: lsgn
        integer :: lt, lix, lcol
        theta(1:nrff) = phase(1:nrff)
        call dgemv('T', norb, nrff, one, omega, norb, Qv, incx, one, theta, incy)
        cosv(:) = cos(theta(:))
        tmp_real(loff+1:loff+nrff) = tmp_real(loff+1:loff+nrff) + lsgn*cosv(1:nrff)*cprod
        if (.not. desc_forces_local) return
        sinv(:) = sin(theta(:))
        call dgemm('T', 'N', norb, 3*m, dimx, one, dQt, dimx, DQmat, MAXX, zero, Bmat, norb)
        call dgemm('T', 'N', nrff, 3*m, norb, one, omega, norb, Bmat, norb, zero, Wmat, nrff)
        do lt = 1, m
          do lix = 1, 3
            lcol = (lt-1)*3 + lix
            gF(1:nrff) = -sinv(1:nrff)*Wmat(1:nrff, lcol)*cprod &
                       + cosv(1:nrff)*cex(lt)*dfcp(lt)*tmp_dxp(lix, jp(lt)+1)/rp(lt)
            tmp_dreal(loff+1:loff+nrff, jp(lt), lix) = tmp_dreal(loff+1:loff+nrff, jp(lt), lix) + lsgn*gF(1:nrff)
          end do
        end do
      end subroutine spip_channel_accum

      ! mean / covariance accumulation of one RFF input sample
      subroutine spip_stats(qs)
        real(kind_double), intent(in) :: qs(:)
        integer :: li, lj
        cnt_n = cnt_n + 1
        if (ltmp_mean_ftnbody) then
          mean_n(1:norb) = mean_n(1:norb) + qs(1:norb)
        else
          do li = 1, norb
            do lj = 1, norb
              covar_n(li, lj) = covar_n(li, lj) + (qs(li) - mean_n(li))*(qs(lj) - mean_n(lj))
            end do
          end do
        end if
      end subroutine spip_stats

      ! ascending insertion sort of the first n entries
      subroutine isort_m(v, n)
        integer, intent(inout) :: v(:)
        integer, intent(in) :: n
        integer :: li, lj, lk
        do li = 1, n - 1
          do lj = li + 1, n
            if (v(lj) < v(li)) then
              lk = v(li); v(li) = v(lj); v(lj) = lk
            end if
          end do
        end do
      end subroutine isort_m

    end subroutine ftnbody_spip_generic

  end module module_compute_body_order



  module module_compute_ftnbody 
  contains 
  subroutine compute_ftnbody(i_start_at, i_final_at, d_n_neigh, d_kind_neigh, iconf)
    use module_kind_variables, ONLY: kind_double
#ifdef MLD_NDM
    use gen_com_m, only : lperiod
    use gen_com_m_ml, only : imm,at,bg
    use tab_imm_m_ml, only : xp
#else
    use ondm_gen_com_m, only : imm, lperiod
    use ondm_tab_imm_m, only : iwmax2, xp
#endif
    use ml_in_ndm_module, only: imm_neigh, desc_forces, rangml
    use derived_types, only: config_real, config_desc
    use module_neigh_local, only: r_cut, & 
                                  max_neigh_local, i_central, i_type, i_type_db, r_central, iw2, &
                                  preallocate_neigh_ja, build_local_neighbours_ja, reallocate_neigh_ja, & 
                                  compute_ur_transformed_distances
    use time_check_general, only: debug_time, MY_MPI_WTIME
    use module_ftnbody, only: tnn_ftbd, t2b_ftbd, t3b_ftbd, t4b_ftbd, t5b_ftbd, tmp_real, tmp_dreal, init_mode_ftnbody, &
                              ftnbody_model_id, FTNBODY_MODEL_GRAMM, FTNBODY_MODEL_CPIP, FTNBODY_MODEL_SPIP
    use module_body_desc, only: l_body_order, dim_desc_body
    use module_compute_body_order, only: ftnbody_order_2, ftnbody_order_3, ftnbody_order_4, ftnbody_order_5, &
                                         ftnbody_gram_order, ftnbody_cpip_order, ftnbody_spip_order
#ifdef MLD_NDM
    use notperiod_mod
#else
    use ondm_transform_coord, only: ondm_notperiod
#endif

    use mld_logger
    use ondm_transform_coord, only: ondm_notperiod
    
    implicit none 
    integer, intent(in)  :: i_start_at, i_final_at
    integer, dimension(imm), intent(out)   :: d_n_neigh
    integer, dimension(imm, imm_neigh), intent(out)    :: d_kind_neigh
    integer, optional    :: iconf

    integer  :: max_neigh
    logical :: desc_forces_local, small 
    real(kind_double), dimension(:, :), allocatable   :: xpnp, tmp_dxp, tmp_xp
    integer :: ja, ja_atom, type_db_ja, icnt, ii, ix, ia  
    real(kind_double) :: t00, t11, t22, t33, t44, t55
    ! real(kind_double), dimension(:), allocatable :: fcut_all, dfcut_all 
    ! real(kind_double) :: fcut, dfcut
    real(kind_double), dimension(:), allocatable  :: ur_central 
    real(kind_double), dimension(:,:), allocatable :: d_ur_central, cos_dxp   
  
    _NAMECURRENT_("compute_ftnbody")

    desc_forces_local = desc_forces .and. (config_real(iconf)%has_force .or. config_real(iconf)%has_stress)
    if (init_mode_ftnbody) desc_forces_local = .false.
    max_neigh = imm_neigh

    if ((i_start_at == 0) .and. (i_final_at == 0)) then
      d_n_neigh(:) = 0
      d_kind_neigh(:, :) = 0
      config_desc(iconf)%energy(:, :) = 0.d0
      return
    end if

    small = .false.
    if (present(iconf)) then
      small = config_real(iconf)%small
    end if
  
    allocate (xpnp(3, imm))
    if (lperiod) then
      xpnp(:, :) = xp(:, :)
    else
#ifdef MLD_NDM
         call notperiod(imm,xp, xpnp,at,bg,.false.)
#else
      call ondm_notperiod(xp, xpnp)
#endif
    end if
    ! call cryst_to_cart (imm, xpnp, bg, -1)
  
    d_n_neigh(:) = 0
    d_kind_neigh(:, :) = 0
    config_desc(iconf)%energy(:, :) = 0.d0
    if (desc_forces_local) config_desc(iconf)%force(:, :, :, :) = 0.d0

    call preallocate_neigh_ja(r_central, i_type, i_central, tmp_dxp, tmp_xp, imm_neigh)

#ifdef MLD_NDM
!JPC    if (i_start_at == 1)  iw2 = 0
!JPC    if (i_start_at >  1)  iw2 = iwmax2(i_start_at - 1)
    iw2=0
#else
    if (i_start_at == 1)  iw2 = 0
    if (i_start_at >  1)  iw2 = iwmax2(i_start_at - 1)  
#endif
    do ja = i_start_at, i_final_at

      ja_atom = ja
      type_db_ja = config_real(iconf)%itype_db(ja)
      if (debug_time) t00 = MY_MPI_WTIME()
      call build_local_neighbours_ja( iconf, ja, imm, xpnp, r_cut, d_n_neigh, d_kind_neigh, &
                                      r_central, i_type, i_type_db, i_central, tmp_dxp, tmp_xp, &
                                      max_neigh_local, iw2)

      if (max_neigh_local == 0) cycle
      if (max_neigh_local > max_neigh) then
        if (rangml == 0) then
          call log_warning("ML: Fatal error in compute_cluster_bonds_and_angles the number of atoms. ")
          call log_warning("ML: The number of max_neig_local is bigger than max_neigh:  r_cut "//vtoa(r_cut)//&
                           " max_neigh_local "// vtoa(max_neigh)//" max_neigh_local "//vtoa(max_neigh_local))
          call log_warning("ML:possible solutions: decrease r_cut or increase max_neigh")
        end if
        call log_critical("ML: decrease rcut for that descriptor")
      end if

      call reallocate_neigh_ja (max_neigh_local, ja_atom, i_central, i_type, i_type_db, r_central, tmp_dxp)

      ! if (allocated(ur_central)) deallocate (ur_central); allocate (ur_central(max_neigh_local))
      ! if (desc_forces_local) then
      !   if (allocated(d_ur_central)) deallocate (d_ur_central); allocate (d_ur_central(3, max_neigh_local))
      ! end if

      ! tmp_dxp is 0-based (central at 0); pass the neighbour section so it aligns
      ! 1-based with r_central inside the transform (d_ur_central(:,k) <-> neighbour k).
      call compute_ur_transformed_distances(desc_forces_local, max_neigh_local, r_central, tmp_dxp(:, 1:max_neigh_local), cos_dxp, ur_central, d_ur_central)
      
      ! do ii = 1, max_neigh_local
      !   call fcut_rij_inout(r_central(ii), r_cut, r_cut_width, r_cut_in, r_cut_width_in, desc_forces_local, fcut, dfcut)
      !   fcut_all(ii) = fcut
      !   dfcut_all(ii) = dfcut
      ! end do

      if (debug_time) then 
        t11 = MY_MPI_WTIME()
        tnn_ftbd = tnn_ftbd + t11 - t00    
      end if
      
      icnt = 0 
      if (l_body_order(2)) then 
        if (allocated(tmp_real)) deallocate (tmp_real); allocate (tmp_real(dim_desc_body(2)))
        if (desc_forces_local) then
          if (allocated(tmp_dreal)) deallocate (tmp_dreal); allocate (tmp_dreal(dim_desc_body(2), max_neigh_local, 3))
        end if 

        select case (ftnbody_model_id)
        case (FTNBODY_MODEL_GRAMM)
          call ftnbody_gram_order(2, type_db_ja, max_neigh_local, i_type_db, r_central, tmp_dxp, &
                                  ur_central, d_ur_central, desc_forces_local)
        case (FTNBODY_MODEL_CPIP)
          call ftnbody_cpip_order(2, type_db_ja, max_neigh_local, i_type_db, r_central, tmp_dxp, &
                                  ur_central, d_ur_central, desc_forces_local)
        case (FTNBODY_MODEL_SPIP)
          call ftnbody_spip_order(2, type_db_ja, max_neigh_local, i_type_db, r_central, tmp_dxp, &
                                  ur_central, d_ur_central, desc_forces_local)
        case default
          call ftnbody_order_2(type_db_ja, max_neigh_local, i_type_db, i_central, r_central, tmp_dxp,   desc_forces_local)
        end select

        ! unit-stride copy (component index innermost on both sides)
        config_desc(iconf)%energy(icnt+1:icnt+dim_desc_body(2), ja) = tmp_real(1:dim_desc_body(2))
        if (desc_forces_local) then
          do ix = 1, 3
            do ia = 1, max_neigh_local
              config_desc(iconf)%force(icnt+1:icnt+dim_desc_body(2), ja, ia, ix) = tmp_dreal(1:dim_desc_body(2), ia, ix)
            end do
          end do
        end if
        icnt = icnt + dim_desc_body(2)
      end if

      if (debug_time) then
        t22 = MY_MPI_WTIME()
        t2b_ftbd = t2b_ftbd + t22 - t11
      end if

      if (l_body_order(3)) then 
        if (allocated(tmp_real)) deallocate (tmp_real); allocate (tmp_real(dim_desc_body(3)))
        if (desc_forces_local) then
          if (allocated(tmp_dreal)) deallocate (tmp_dreal); allocate (tmp_dreal(dim_desc_body(3), max_neigh_local, 3))
        end if 

        select case (ftnbody_model_id)
        case (FTNBODY_MODEL_GRAMM)
          call ftnbody_gram_order(3, type_db_ja, max_neigh_local, i_type_db, r_central, tmp_dxp, &
                                  ur_central, d_ur_central, desc_forces_local)
        case (FTNBODY_MODEL_CPIP)
          call ftnbody_cpip_order(3, type_db_ja, max_neigh_local, i_type_db, r_central, tmp_dxp, &
                                  ur_central, d_ur_central, desc_forces_local)
        case (FTNBODY_MODEL_SPIP)
          call ftnbody_spip_order(3, type_db_ja, max_neigh_local, i_type_db, r_central, tmp_dxp, &
                                  ur_central, d_ur_central, desc_forces_local)
        case default
          call ftnbody_order_3(type_db_ja, max_neigh_local, i_type_db, i_central, r_central, tmp_dxp, ur_central, d_ur_central, desc_forces_local)
        end select
        ! unit-stride copy (component index innermost on both sides)
        config_desc(iconf)%energy(icnt+1:icnt+dim_desc_body(3), ja) = tmp_real(1:dim_desc_body(3))
        if (desc_forces_local) then
          do ix = 1, 3
            do ia = 1, max_neigh_local
              config_desc(iconf)%force(icnt+1:icnt+dim_desc_body(3), ja, ia, ix) = tmp_dreal(1:dim_desc_body(3), ia, ix)
            end do
          end do
        end if
        icnt = icnt + dim_desc_body(3)

      end if

      if (debug_time) then 
        t33 = MY_MPI_WTIME()
        t3b_ftbd = t3b_ftbd + t33 - t22    
      end if

      if (l_body_order(4)) then 
        if (allocated(tmp_real)) deallocate (tmp_real); allocate (tmp_real(dim_desc_body(4)))
        if (desc_forces_local) then
          if (allocated(tmp_dreal)) deallocate (tmp_dreal); allocate (tmp_dreal(dim_desc_body(4), max_neigh_local, 3))
        end if 

        select case (ftnbody_model_id)
        case (FTNBODY_MODEL_GRAMM)
          call ftnbody_gram_order(4, type_db_ja, max_neigh_local, i_type_db, r_central, tmp_dxp, &
                                  ur_central, d_ur_central, desc_forces_local)
        case (FTNBODY_MODEL_CPIP)
          call ftnbody_cpip_order(4, type_db_ja, max_neigh_local, i_type_db, r_central, tmp_dxp, &
                                  ur_central, d_ur_central, desc_forces_local)
        case (FTNBODY_MODEL_SPIP)
          call ftnbody_spip_order(4, type_db_ja, max_neigh_local, i_type_db, r_central, tmp_dxp, &
                                  ur_central, d_ur_central, desc_forces_local)
        case default
          call ftnbody_order_4(type_db_ja, max_neigh_local, i_type_db, i_central, r_central, tmp_dxp, ur_central, d_ur_central, desc_forces_local)
        end select
        ! unit-stride copy (component index innermost on both sides)
        config_desc(iconf)%energy(icnt+1:icnt+dim_desc_body(4), ja) = tmp_real(1:dim_desc_body(4))
        if (desc_forces_local) then
          do ix = 1, 3
            do ia = 1, max_neigh_local
              config_desc(iconf)%force(icnt+1:icnt+dim_desc_body(4), ja, ia, ix) = tmp_dreal(1:dim_desc_body(4), ia, ix)
            end do
          end do
        end if
        icnt = icnt + dim_desc_body(4)
      end if


      if (debug_time) then 
        t44 = MY_MPI_WTIME()
        t4b_ftbd = t4b_ftbd + t44 - t33    
      end if

      if (l_body_order(5)) then 
        if (allocated(tmp_real)) deallocate (tmp_real); allocate (tmp_real(dim_desc_body(5)))
        if (desc_forces_local) then
          if (allocated(tmp_dreal)) deallocate (tmp_dreal); allocate (tmp_dreal(dim_desc_body(5), max_neigh_local, 3))
        end if 

        select case (ftnbody_model_id)
        case (FTNBODY_MODEL_GRAMM)
          call ftnbody_gram_order(5, type_db_ja, max_neigh_local, i_type_db, r_central, tmp_dxp, &
                                  ur_central, d_ur_central, desc_forces_local)
        case (FTNBODY_MODEL_CPIP)
          call ftnbody_cpip_order(5, type_db_ja, max_neigh_local, i_type_db, r_central, tmp_dxp, &
                                  ur_central, d_ur_central, desc_forces_local)
        case (FTNBODY_MODEL_SPIP)
          call ftnbody_spip_order(5, type_db_ja, max_neigh_local, i_type_db, r_central, tmp_dxp, &
                                  ur_central, d_ur_central, desc_forces_local)
        case default
          call ftnbody_order_5(type_db_ja, max_neigh_local, i_type_db, i_central, r_central, tmp_dxp, ur_central, d_ur_central,  desc_forces_local)
        end select
        ! unit-stride copy (component index innermost on both sides)
        config_desc(iconf)%energy(icnt+1:icnt+dim_desc_body(5), ja) = tmp_real(1:dim_desc_body(5))
        if (desc_forces_local) then
          do ix = 1, 3
            do ia = 1, max_neigh_local
              config_desc(iconf)%force(icnt+1:icnt+dim_desc_body(5), ja, ia, ix) = tmp_dreal(1:dim_desc_body(5), ia, ix)
            end do
          end do
        end if
        icnt = icnt + dim_desc_body(5)
      end if

      if (debug_time) then 
        t55 = MY_MPI_WTIME()
        t5b_ftbd = t5b_ftbd + t55 - t44    
      end if

      if (desc_forces_local) then
        do ia = 1, max_neigh_local
          config_desc(iconf)%force(:, ja, 0, :) = config_desc(iconf)%force(:, ja, 0, :) - config_desc(iconf)%force(:, ja, ia, :)
        end do
      end if
      
    end do ! end_ja   

    _MLD_BEGIN_

    _MLD_END_ 

  end subroutine compute_ftnbody
  end module module_compute_ftnbody 





