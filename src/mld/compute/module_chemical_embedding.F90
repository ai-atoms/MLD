!> @brief Module for chemical embedding in CHEMMAP_HSVD
!> @details Implements PCA-based chemical embeddings from physical properties
!>          Following section 1.6 of MiladyNoteTechnique4.pdf
!> @author MCM
!> @date December 2025
module module_chemical_embedding
  use iso_fortran_env, only: dp => real64
  use mld_logger
  use mld_string, only: vtoa
  implicit none
  private

  !> Chemical properties for each element
  type :: type_chemical_properties
    real(dp) :: atomic_radius      ! Covalent radius (Å)
    real(dp) :: electronegativity  ! Pauling electronegativity
    real(dp) :: cohesive_energy    ! Cohesive energy (eV)
    real(dp) :: valence_electrons  ! Number of valence electrons
    real(dp) :: atomic_mass        ! Atomic mass (amu)
    real(dp) :: density            ! Elemental density (g/cm³)
  end type type_chemical_properties

  !> Chemical embedding container
  type :: type_chemical_embedding
    integer :: n_species                              ! Number of species
    integer :: d_prop                                 ! Number of properties (6)
    integer :: d_e                                    ! Embedding dimension
    real(dp), dimension(:,:), allocatable :: U_properties  ! (n_species, d_prop)
    real(dp), dimension(:), allocatable :: u_mean     ! Mean properties (d_prop)
    real(dp), dimension(:,:), allocatable :: V_basis  ! PCA basis (d_prop, d_prop)
    real(dp), dimension(:), allocatable :: sigma      ! Singular values
    real(dp), dimension(:,:), allocatable :: embeddings ! (n_species, d_e)
    
    ! Isotropic and anisotropic maps
    real(dp), dimension(:,:), allocatable :: W_iso        ! (d_p, 2*d_e)
    real(dp), dimension(:,:), allocatable :: W_aniso      ! (d_p, 2*d_e)
    real(dp), dimension(:,:,:), allocatable :: W_l        ! (lmax, d_p, 2*d_e)
    
    ! Parameters
    real(dp) :: beta_gamma                            ! Decay parameter for gamma_l
    integer :: d_p                                    ! Pair descriptor dimension
    integer :: n_aniso_modes                          ! Number of anisotropic modes
    
    logical :: is_initialized = .false.
  contains
    procedure :: init => init_chemical_embedding
    procedure :: build_embeddings => build_chemical_embeddings
    procedure :: build_isotropic_map => build_W_iso
    procedure :: build_anisotropic_map => build_W_aniso
    procedure :: build_l_dependent_maps => build_W_l
    procedure :: get_pair_descriptor => get_c_pair
    procedure :: cleanup => cleanup_chemical_embedding
  end type type_chemical_embedding

  ! Global instance
  type(type_chemical_embedding), save :: chem_embed

  public :: type_chemical_embedding, chem_embed
  public :: init_chemmap_properties, get_chemical_pair_coeff, get_c_pair_wrapper
  public :: write_chemical_embedding_debug

contains

  !> Initialize chemical properties for common elements
  subroutine init_chemmap_properties(species_list, n_species)
    use module_chemical_species, only: fix_ch_elements
    integer, intent(in) :: n_species
    character(len=2), dimension(:), intent(in) :: species_list
    type(type_chemical_properties) :: props
    integer :: i, j
    
    integer :: d_e, d_p
    real(dp) :: d_p_factor
    
    call log_info("ML: Initializing chemical embedding for CHEMMAP_HSVD")
    call log_info("ML: Number of species: "//vtoa(n_species))
    
    ! Initialize the embedding structure
    ! d_p_factor controls projection dimension relative to embedding dimension
    ! Options: 1.0 (minimal), 1.8 (intermediate), 2.0 (full sum+contrast)
    d_e = 3
    d_p_factor = 1.8_dp
    d_p = nint(d_p_factor * real(d_e, dp))
    
    
    call log_info("ML: d_e = "//vtoa(d_e)//" (embedding dimension)")
    call log_info("ML: d_p = "//vtoa(d_p)//" (projection dimension, factor="//vtoa(d_p_factor,'f4.2')//")")
    
    ! Set n_aniso to use the first few principal chemical directions (typically 1-3)
    ! This determines how many SVD modes are used for chemical contrast
    !call chem_embed%init(n_species, d_e=d_e, d_p=d_p, n_aniso=min(d_p, d_e), beta=0.5_dp)
    call chem_embed%init(n_species, d_e=d_e, d_p=d_p, n_aniso=2, beta=0.5_dp)
    
    ! Fill properties for each species
    do i = 1, n_species
      props = get_element_properties(trim(species_list(i)))
      chem_embed%U_properties(i, 1) = props%atomic_radius
      chem_embed%U_properties(i, 2) = props%electronegativity
      chem_embed%U_properties(i, 3) = props%cohesive_energy
      chem_embed%U_properties(i, 4) = props%valence_electrons
      chem_embed%U_properties(i, 5) = props%atomic_mass
      chem_embed%U_properties(i, 6) = props%density
      
      call log_info("ML:   "//trim(species_list(i))//": r="//vtoa(props%atomic_radius,'f7.3')// &
                   " chi="//vtoa(props%electronegativity,'f6.3')//" E_coh="//vtoa(props%cohesive_energy,'f8.3'))
    end do
    
    ! Build embeddings via PCA
    call chem_embed%build_embeddings()
    
    ! Build isotropic and anisotropic maps
    call chem_embed%build_isotropic_map()
    call chem_embed%build_anisotropic_map()
    call chem_embed%build_l_dependent_maps(lmax=6)
    
    call log_info("ML: Chemical embedding initialization complete")
  end subroutine init_chemmap_properties

  !> Initialize the chemical embedding structure
  subroutine init_chemical_embedding(this, n_species, d_e, d_p, n_aniso, beta)
    class(type_chemical_embedding), intent(inout) :: this
    integer, intent(in) :: n_species
    integer, intent(in) :: d_e          ! Embedding dimension
    integer, intent(in) :: d_p          ! Pair descriptor dimension
    integer, intent(in) :: n_aniso      ! Number of anisotropic modes
    real(dp), intent(in) :: beta        ! Decay parameter
    
    this%n_species = n_species
    this%d_prop = 6  ! Fixed: radius, chi, E_coh, Z, mass, density
    this%d_e = min(d_e, this%d_prop)  ! Cannot exceed number of properties
    this%d_p = d_p
    this%n_aniso_modes = n_aniso
    this%beta_gamma = beta
    
    ! Allocate arrays
    allocate(this%U_properties(n_species, this%d_prop))
    allocate(this%u_mean(this%d_prop))
    allocate(this%V_basis(this%d_prop, this%d_prop))
    allocate(this%sigma(this%d_prop))
    allocate(this%embeddings(n_species, this%d_e))
    allocate(this%W_iso(d_p, 2*this%d_e))
    allocate(this%W_aniso(d_p, 2*this%d_e))
    
    this%U_properties = 0.0_dp
    this%u_mean = 0.0_dp
    this%V_basis = 0.0_dp
    this%sigma = 0.0_dp
    this%embeddings = 0.0_dp
    this%W_iso = 0.0_dp
    this%W_aniso = 0.0_dp
    
    this%is_initialized = .true.
  end subroutine init_chemical_embedding

  !> Build chemical embeddings via PCA/SVD
  subroutine build_chemical_embeddings(this)
    class(type_chemical_embedding), intent(inout) :: this
    real(dp), dimension(:,:), allocatable :: U_centered, U_copy
    real(dp), dimension(:,:), allocatable :: U_svd, Vt_svd
    real(dp), dimension(:), allocatable :: u_std
    integer :: i, j, info
    
    ! Compute mean
    do j = 1, this%d_prop
      this%u_mean(j) = sum(this%U_properties(:, j)) / real(this%n_species, dp)
    end do
    
    ! Compute standard deviation
    allocate(u_std(this%d_prop))
    do j = 1, this%d_prop
      u_std(j) = sqrt(sum((this%U_properties(:, j) - this%u_mean(j))**2) / real(this%n_species, dp))
      if (u_std(j) < 1.0e-12_dp) u_std(j) = 1.0_dp  ! Avoid division by zero
    end do
    
    ! Center and standardize the data
    allocate(U_centered(this%n_species, this%d_prop))
    do i = 1, this%n_species
      do j = 1, this%d_prop
        U_centered(i, j) = (this%U_properties(i, j) - this%u_mean(j)) / u_std(j)
      end do
    end do
    
    ! Perform SVD: U_centered = U * Sigma * V^T
    ! For embedding we need V (right singular vectors)
    allocate(U_copy(this%n_species, this%d_prop))
    allocate(U_svd(this%n_species, this%n_species))
    allocate(Vt_svd(this%d_prop, this%d_prop))
    U_copy = U_centered
    
    ! Call LAPACK SVD
    call compute_svd_wrapper(U_copy, U_svd, this%sigma, Vt_svd, info)
    
    if (info /= 0) then
      call log_error("ML: SVD failed in chemical embedding with info="//vtoa(info))
      stop
    end if
    
    ! V_basis is the transpose of VT (eigvecVT is already V^T from SVD)
    do i = 1, this%d_prop
      do j = 1, this%d_prop
        this%V_basis(i, j) = svd_solver%eigvecVT(j, i)
      end do
    end do
    
    call svd_solver%destroy()
    
    ! Compute embeddings: e(mu) = V^T * (u(mu) - u_mean)
    ! embeddings(mu, :) = V_basis^T(:, 1:d_e) * U_centered(mu, :)
    do i = 1, this%n_species
      do j = 1, this%d_e
        this%embeddings(i, j) = dot_product(this%V_basis(:, j), U_centered(i, :))
      end do
    end do
    
    call log_info("ML: Chemical embedding: sigma(1:3)="//vtoa(this%sigma(1),'f10.4')//" "// &
                 vtoa(this%sigma(2),'f10.4')//" "//vtoa(this%sigma(3),'f10.4'))
    
    ! Check rank of chemical property matrix
    if (this%d_e > this%n_species - 1) then
      call log_warning("ML: d_e="//vtoa(this%d_e)//" exceeds max rank="//vtoa(this%n_species-1)// &
                      " for "//vtoa(this%n_species)//" species")
      call log_warning("ML: Effective chemical embedding dimension limited to "//vtoa(this%n_species-1))
    end if
    
    deallocate(U_centered, u_std)
  end subroutine build_chemical_embeddings

  !> Build isotropic map W_0 (pair-sum)
  subroutine build_W_iso(this)
    class(type_chemical_embedding), intent(inout) :: this
    integer :: p, e
    

    this%W_iso = 0.0_dp
    do p = 1, min(this%d_p, this%d_e) 
      this%W_iso(p, p) = 1.0_dp              
      this%W_iso(p, this%d_e + p) = 1.0_dp  
    end do
    
    call log_info("ML: Built isotropic map W_0 (pair-sum) with d_p="//vtoa(this%d_p)//" d_e="//vtoa(this%d_e))
  end subroutine build_W_iso

  !> Build anisotropic map W_aniso (pair-contrast)
  !> Uses eigenvectors from PCA to create contrast projections
  !> Row structure: [v_k, -v_k] where v_k is the k-th principal direction
  subroutine build_W_aniso(this)
    class(type_chemical_embedding), intent(inout) :: this
    integer :: p, mode, e
    
    ! W_aniso structure: each row projects onto differences along principal directions
    ! Row p uses mode = mod(p-1, n_aniso_modes) + 1
    ! The contrast is built from the eigenvector components: [v_mode, -v_mode]
    this%W_aniso = 0.0_dp
    do p = 1, this%d_p
      mode = mod(p-1, min(this%n_aniso_modes, this%d_e)) + 1
      
      ! For this mode, use the corresponding eigenvector components
      ! First d_e columns: +v_mode (for mu_a embedding)
      ! Last d_e columns: -v_mode (for mu_j embedding) - creates contrast
      do e = 1, this%d_e
        this%W_aniso(p, e) = this%V_basis(e, mode)              ! +v_mode component for mu_a
        this%W_aniso(p, this%d_e + e) = -this%V_basis(e, mode)  ! -v_mode component for mu_j
      end do
    end do
    
    call log_info("ML: Built anisotropic map W_aniso with d_p="//vtoa(this%d_p)//" using "//vtoa(this%n_aniso_modes)//" modes")
  end subroutine build_W_aniso

  !> Build l-dependent maps: W^(l) = gamma_l * W_0 + (1 - gamma_l) * W_aniso
  subroutine build_W_l(this, lmax)
    class(type_chemical_embedding), intent(inout) :: this
    integer, intent(in) :: lmax
    real(dp) :: gamma_l
    integer :: l, i, j
    
    if (allocated(this%W_l)) deallocate(this%W_l)
    allocate(this%W_l(0:lmax, this%d_p, 2*this%d_e))
    
    do l = 0, lmax
      gamma_l = exp(-this%beta_gamma * real(l, dp))
      
      do i = 1, this%d_p
        do j = 1, 2*this%d_e
          this%W_l(l, i, j) = gamma_l * this%W_iso(i, j) + (1.0_dp - gamma_l) * this%W_aniso(i, j)
        end do
      end do
    end do
    
    call log_info("ML: Built l-dependent maps for l=0 to "//vtoa(lmax))
    call log_info("ML:   gamma(l=0)="//vtoa(exp(0.0_dp),'f8.4')//" gamma(l=2)="//vtoa(exp(-2*this%beta_gamma),'f8.4'))
  end subroutine build_W_l

  !> Get pair chemical descriptor c^(l)_p(mu_a, mu_j)
  function get_c_pair(this, mu_a, mu_j, l) result(c_pair)
    class(type_chemical_embedding), intent(in) :: this
    integer, intent(in) :: mu_a, mu_j, l
    real(dp), dimension(this%d_p) :: c_pair
    real(dp), dimension(2*this%d_e) :: z_pair
    integer :: i, j
    
    ! Concatenate embeddings: z = [e(mu_a), e(mu_j)]
    z_pair(1:this%d_e) = this%embeddings(mu_a, :)
    z_pair(this%d_e+1:2*this%d_e) = this%embeddings(mu_j, :)
    
    ! c^(l) = W^(l) * z
    c_pair = 0.0_dp
    do i = 1, this%d_p
      do j = 1, 2*this%d_e
        c_pair(i) = c_pair(i) + this%W_l(l, i, j) * z_pair(j)
      end do
    end do
  end function get_c_pair

  !> Wrapper to get single chemical coefficient for specific mu_j index
  function get_c_pair_wrapper(mu_a, mu_j, l, p_idx) result(c_val)
    integer, intent(in) :: mu_a, mu_j, l, p_idx
    real(dp) :: c_val
    real(dp), dimension(chem_embed%d_p) :: c_pair
    
    if (.not. chem_embed%is_initialized) then
      c_val = 0.0_dp
      if (mu_a == mu_j) c_val = 1.0_dp  ! Fallback to delta
      return
    end if
    
    c_pair = get_c_pair(chem_embed, mu_a, mu_j, l)
    
    if (p_idx >= 1 .and. p_idx <= chem_embed%d_p) then
      c_val = c_pair(p_idx)
    else
      c_val = 0.0_dp
    end if
  end function get_c_pair_wrapper

  !> Public interface to get chemical coefficient for pair
  subroutine get_chemical_pair_coeff(mu_a, mu_j, l, n, c_np)
    integer, intent(in) :: mu_a, mu_j, l, n
    real(dp), intent(out) :: c_np
    real(dp), dimension(chem_embed%d_p) :: c_pair
    
    if (.not. chem_embed%is_initialized) then
      c_np = 1.0_dp  ! Fallback
      return
    end if
    
    c_pair = get_c_pair(chem_embed, mu_a, mu_j, l)
    
    ! Map to specific (n,p) index
    ! For now, simple mapping: p = mod(n-1, d_p) + 1
    c_np = c_pair(mod(n-1, chem_embed%d_p) + 1)
  end subroutine get_chemical_pair_coeff

  !> Get element properties from database
  function get_element_properties(element_symbol) result(props)
    character(len=*), intent(in) :: element_symbol
    type(type_chemical_properties) :: props
    
    ! Database of common elements
    ! Format: radius(Å), electronegativity, cohesive_energy(eV), valence_e, mass(amu), density(g/cm³)
    
    select case(trim(element_symbol))
      case('H')
        props = type_chemical_properties(0.31_dp, 2.20_dp, 0.00_dp, 1.0_dp, 1.008_dp, 0.09_dp)
      case('C')
        props = type_chemical_properties(0.76_dp, 2.55_dp, 7.37_dp, 4.0_dp, 12.011_dp, 2.27_dp)
      case('N')
        props = type_chemical_properties(0.71_dp, 3.04_dp, 0.00_dp, 5.0_dp, 14.007_dp, 1.25_dp)
      case('O')
        props = type_chemical_properties(0.66_dp, 3.44_dp, 0.00_dp, 6.0_dp, 15.999_dp, 1.43_dp)
      case('Fe')
        props = type_chemical_properties(1.32_dp, 1.83_dp, 4.28_dp, 8.0_dp, 55.845_dp, 7.87_dp)
      case('Ni')
        props = type_chemical_properties(1.24_dp, 1.91_dp, 4.44_dp, 10.0_dp, 58.693_dp, 8.91_dp)
      case('Cu')
        props = type_chemical_properties(1.32_dp, 1.90_dp, 3.49_dp, 11.0_dp, 63.546_dp, 8.96_dp)
      case('Ti')
        props = type_chemical_properties(1.60_dp, 1.54_dp, 4.85_dp, 4.0_dp, 47.867_dp, 4.51_dp)
      case('Ta')
        props = type_chemical_properties(1.46_dp, 1.50_dp, 8.10_dp, 5.0_dp, 180.948_dp, 16.65_dp)
      case('W')
        props = type_chemical_properties(1.41_dp, 2.36_dp, 8.90_dp, 6.0_dp, 183.840_dp, 19.25_dp)
      case('Nb')
        props = type_chemical_properties(1.64_dp, 1.60_dp, 7.57_dp, 5.0_dp, 92.906_dp, 8.57_dp)
      case('Mo')
        props = type_chemical_properties(1.54_dp, 2.16_dp, 6.82_dp, 6.0_dp, 95.950_dp, 10.28_dp)
      case('V')
        props = type_chemical_properties(1.53_dp, 1.63_dp, 5.31_dp, 5.0_dp, 50.942_dp, 6.11_dp)
      case('Cr')
        props = type_chemical_properties(1.39_dp, 1.66_dp, 4.10_dp, 6.0_dp, 51.996_dp, 7.19_dp)
      case('Al')
        props = type_chemical_properties(1.21_dp, 1.61_dp, 3.39_dp, 3.0_dp, 26.982_dp, 2.70_dp)
      case('Si')
        props = type_chemical_properties(1.11_dp, 1.90_dp, 4.63_dp, 4.0_dp, 28.085_dp, 2.33_dp)
      case default
        ! Default values for unknown elements
        call log_warning("ML: Unknown element '"//trim(element_symbol)//"', using default properties")
        props = type_chemical_properties(1.5_dp, 2.0_dp, 5.0_dp, 5.0_dp, 100.0_dp, 5.0_dp)
    end select
  end function get_element_properties

  !> Cleanup
  subroutine cleanup_chemical_embedding(this)
    class(type_chemical_embedding), intent(inout) :: this
    
    if (allocated(this%U_properties)) deallocate(this%U_properties)
    if (allocated(this%u_mean)) deallocate(this%u_mean)
    if (allocated(this%V_basis)) deallocate(this%V_basis)
    if (allocated(this%sigma)) deallocate(this%sigma)
    if (allocated(this%embeddings)) deallocate(this%embeddings)
    if (allocated(this%W_iso)) deallocate(this%W_iso)
    if (allocated(this%W_aniso)) deallocate(this%W_aniso)
    if (allocated(this%W_l)) deallocate(this%W_l)
    
    this%is_initialized = .false.
  end subroutine cleanup_chemical_embedding

  !> Write chemical embedding debug information to files
  subroutine write_chemical_embedding_debug(species_list, lmax)
    use mld_mpi, only: mld_rank
    character(len=2), dimension(:), intent(in) :: species_list
    integer, intent(in) :: lmax
    integer :: mu_a, mu_j, l, p, i, unit_emb, unit_coeff
    real(dp), dimension(:), allocatable :: c_pair
    character(len=256) :: filename
    
    if (mld_rank /= 0) return  ! Only rank 0 writes
    if (.not. chem_embed%is_initialized) then
      call log_warning("ML: Chemical embedding not initialized, skipping debug output")
      return
    end if
    
    ! Write embeddings e(μ) for each species
    filename = 'chemmap_embeddings.dat'
    open(newunit=unit_emb, file=trim(filename), status='replace', action='write')
    write(unit_emb, '(A)') '# Chemical embeddings e(mu) for each species'
    write(unit_emb, '(A,I0,A,I0)') '# n_species=', chem_embed%n_species, ' d_e=', chem_embed%d_e
    write(unit_emb, '(A)') '# Format: mu species e_1 e_2 ... e_d_e'
    do mu_a = 1, chem_embed%n_species
      write(unit_emb, '(I3,1X,A2,20(1X,ES14.6))') mu_a, species_list(mu_a), &
                                                   (chem_embed%embeddings(mu_a, i), i=1, chem_embed%d_e)
    end do
    close(unit_emb)
    call log_info("ML: Wrote chemical embeddings to "//trim(filename))
    
    ! Write chemical coefficients c^(l)_p(μ_a, μ_j) for all pairs and l
    filename = 'chemmap_coefficients.dat'
    open(newunit=unit_coeff, file=trim(filename), status='replace', action='write')
    write(unit_coeff, '(A)') '# Chemical coefficients c^(l)_p(mu_a, mu_j)'
    write(unit_coeff, '(A,I0,A,I0,A,I0)') '# n_species=', chem_embed%n_species, &
                                          ' d_p=', chem_embed%d_p, ' lmax=', lmax
    write(unit_coeff, '(A)') '# Format: mu_a mu_j l p c_p'
    
    allocate(c_pair(chem_embed%d_p))
    do l = 0, lmax
      do mu_a = 1, chem_embed%n_species
        do mu_j = 1, chem_embed%n_species
          c_pair = get_c_pair(chem_embed, mu_a, mu_j, l)
          do p = 1, chem_embed%d_p
            write(unit_coeff, '(4I4,1X,ES14.6)') mu_a, mu_j, l, p, c_pair(p)
          end do
        end do
      end do
    end do
    deallocate(c_pair)
    close(unit_coeff)
    call log_info("ML: Wrote chemical coefficients to "//trim(filename))
    
    ! Write gamma(l) decay for reference
    filename = 'chemmap_gamma.dat'
    open(newunit=unit_coeff, file=trim(filename), status='replace', action='write')
    write(unit_coeff, '(A)') '# l-dependent mixing parameter gamma(l) = exp(-beta*l)'
    write(unit_coeff, '(A,ES14.6)') '# beta=', chem_embed%beta_gamma
    write(unit_coeff, '(A)') '# Format: l gamma(l)'
    do l = 0, lmax
      write(unit_coeff, '(I4,1X,ES14.6)') l, exp(-chem_embed%beta_gamma * l)
    end do
    close(unit_coeff)
    call log_info("ML: Wrote gamma(l) to "//trim(filename))
    
  end subroutine write_chemical_embedding_debug

end module module_chemical_embedding
