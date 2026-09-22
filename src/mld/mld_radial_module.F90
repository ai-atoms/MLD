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

module RadialFunctions
  use module_kind_variables, only: kind_double
  use ml_in_ndm_module, only: one_pi
  use mld_logger
  implicit none

  private
  public :: RadialFunction

  type RadialFunction
    !private 
    !public :: radial, d_radial 

    logical :: check_init=.false. 
    real(kind_double), allocatable :: values(:)
    integer :: nsp, nradial 
    character(10) :: radial_type


    real(kind_double), dimension(:), allocatable :: radial, d_radial
    real(kind_double), dimension(:), allocatable :: only_radial, d_only_radial
    real(kind_double) :: r_cut_in, r_cut_out, delta_cut_in, delta_cut_out 
    integer :: n_elem 

    ! Chebyshev polynomials as radial parts. 
    real(kind_double), dimension(:), allocatable :: chebT, chebU, d_chebT, d_chebU 
    integer :: n_cheb
    real(kind_double) :: lambda

    ! simplified Bessel polynomials as radial parts.
    integer :: n_bessel
    real(kind_double), dimension(:), allocatable  :: d_bess, e_bess
    real(kind_double), dimension(:), allocatable  :: bf, bg, d_bf, d_bg


  contains
    procedure :: init
    procedure :: evaluate => evaluate_radial
    procedure :: chebTU
    procedure :: radial_powTcheb, radial_expTcheb 
  end type RadialFunction


  contains

  ! Initialize the RadialFunction object with parameters
  !subroutine init(this, radial_type, params)
  subroutine init(this, int_radial_type, r_cut_in, r_cut, r_cut_width_in,  r_cut_width, lambda, acenmax)  
    use mld_mpi, only: mld_mpi_abort
    class(RadialFunction), intent(inout) :: this
    integer, intent(in) :: int_radial_type
    real(kind_double), intent(in) :: r_cut_in, r_cut, r_cut_width_in,  r_cut_width, lambda
    integer, intent(in) :: acenmax 
    integer :: ii 
    real(kind_double) :: rzozo 
    ! function type can have values: 'powTcheb', 'expTcheb', 'simpBessel', 'powWave'
    !powTcheb: pow * Tchebyshev from 0 to n_cheb-1 from R. Drautz et al. PRB 72, 104108 (2005)
    !expTcheb: exp * Tchebyshev from 0 to n_cheb-1 from R. Drautz and Gironcolli et al. 
    !simpBessel comes from Kocer et al. 

    rzozo = r_cut_width
    if (int_radial_type == 1) then
      this%radial_type = 'powTcheb'
    else if (int_radial_type == 2) then
      this%radial_type = 'expTcheb'
    else if (int_radial_type == 3) then
      this%radial_type = 'simpBessel'
    else if (int_radial_type == 4) then
      this%radial_type = 'powWave'
    else
      call log_critical('Unknown radial function type for ACE-like descriptor!')
      stop 'stop in init RadialFunctions'
    end if


    this%r_cut_out = r_cut          ! Cutoff radius for the radial function
    this%r_cut_in = r_cut_in        ! Cutoff radius for the inner part of the radial function
    this%delta_cut_in = r_cut_width_in  ! Width of the inner part of the radial function
    this%delta_cut_out = r_cut_width    ! Width of the outer part of the radial function
    this%n_elem = 1                 ! number of elements in the radial function
    this%lambda = lambda            ! lambda from the radial function

    if (this%radial_type == 'simpBessel') then
      !TODOace 
      this%n_bessel = acenmax 
      ! normaly are n_bessel + 1 because the radial function is defined from 0 to n_bessel
      ! but radial 0 is consrtant value 1, so we need only n_bessel
      !this%nradial = this%n_bessel + 1 
      this%nradial = this%n_bessel  
      if (allocated(this%d_bess)) deallocate(this%d_bess) ; allocate(this%d_bess(0:this%n_bessel))
      if (allocated(this%e_bess)) deallocate(this%e_bess) ; allocate(this%e_bess(0:this%n_bessel))
      this%d_bess(0) = 1.d0
      this%e_bess(0) = 0.d0
      do ii = 1, this%n_bessel 
        this%e_bess(ii) = dble((ii * (ii+2))**2)/dble(4.d0 * (ii+1)**4 + 1)
        this%d_bess(ii) = 1.d0 - this%e_bess(ii) / this%d_bess(ii-1)
      end do

      if (allocated(this%bf)) deallocate(this%bf) ; allocate(this%bf(0:this%n_bessel))
      if (allocated(this%bg)) deallocate(this%bg) ; allocate(this%bg(0:this%n_bessel))
      if (allocated(this%d_bf)) deallocate(this%d_bf) ; allocate(this%d_bf(0:this%n_bessel))
      if (allocated(this%d_bg)) deallocate(this%d_bg) ; allocate(this%d_bg(0:this%n_bessel))
    else if (this%radial_type == 'powTcheb') then
      if (this%lambda <= 1.0 ) then
         call log_critical('Big error in init of radial part. lambda should be bigger than 1 for powTcheb (ace_radial_type=2). Fatal. ')
         call mld_mpi_abort("stop in init RadialFunctions lambda too small.")
      end if
      this%n_cheb = acenmax 
      if (allocated(this%chebT)) deallocate(this%chebT)     ; allocate(this%chebT(0:this%n_cheb))
      if (allocated(this%chebU)) deallocate(this%chebU)     ; allocate(this%chebU(0:this%n_cheb))
      if (allocated(this%d_chebT)) deallocate(this%d_chebT) ; allocate(this%d_chebT(0:this%n_cheb))
      if (allocated(this%d_chebU)) deallocate(this%d_chebU) ; allocate(this%d_chebU(0:this%n_cheb))
      ! normaly us n_cheb + 1 because the radial function is defined from 0 to n_cheb
      ! but radial 0 is constant value 1, so we need only n_cheb
      !this%nradial = this%n_cheb + 1 
      this%nradial = this%n_cheb  
      if (this%lambda <= 1 ) then
         call log_critical('Big error in init of radial part. lambda should be bigger than 1. Fatal. ')
         stop 'stop in init RadialFunctions'
      end if 

    else if (this%radial_type == 'expTcheb') then
      this%n_cheb = acenmax 
      if (allocated(this%chebT)) deallocate(this%chebT)     ;  allocate(this%chebT(0:this%n_cheb))
      if (allocated(this%chebU)) deallocate(this%chebU)     ;  allocate(this%chebU(0:this%n_cheb))
      if (allocated(this%d_chebT)) deallocate(this%d_chebT) ; allocate(this%d_chebT(0:this%n_cheb))
      if (allocated(this%d_chebU)) deallocate(this%d_chebU) ;  allocate(this%d_chebU(0:this%n_cheb))
      ! normaly us n_cheb + 1 because the radial function is defined from 0 to n_cheb
      ! but radial 0 is constant value 1, so we need only n_cheb
      !this%nradial = this%n_cheb + 1 
      this%nradial = this%n_cheb 

    else if (this%radial_type == 'powWave') then
      !TODOace
      this%nradial = this%n_cheb
      call log_critical('powWave not implemented yet in RadialFunctions!')
      stop 'stop in init RadialFunctions wave' 

    else
      call log_critical("Unknown function type!")
      stop
    end if


    if (allocated(this%only_radial)) deallocate(this%only_radial) ; allocate(this%only_radial(0:this%nradial))
    if (allocated(this%d_only_radial)) deallocate(this%d_only_radial) ; allocate(this%d_only_radial(0:this%nradial))

    if (allocated(this%radial)) deallocate(this%radial) ; allocate(this%radial(0:this%nradial))
    if (allocated(this%d_radial)) deallocate(this%d_radial) ; allocate(this%d_radial(0:this%nradial))

    this%check_init = .true. 
  end subroutine init


  subroutine chebTU(this,  xx)
    ! provide the chebyshev polynomials and their derivatives
    implicit none
    class(RadialFunction), intent(inout) :: this
    real(kind_double), intent(in)    :: xx
    !local variables 
    integer  :: nn

    this%chebT(0) = 1.d0
    this%chebT(1) = xx
    this%chebU(0) = 1.d0
    this%chebU(1) = 2.d0*xx

    this%d_chebT(0) = 0.d0
    this%d_chebT(1) = 1.d0
    if (this%n_cheb .ge. 2) then
      do nn = 2, this%n_cheb
        this%chebT(nn) = xx*this%chebT(nn - 1) - (1.d0 - xx**2)*this%chebU(nn - 2)
        this%chebU(nn) = xx*this%chebU(nn - 1) + this%chebT(nn)
      end do

      do nn = 2, this%n_cheb
        this%d_chebT(nn) = dble(nn)*this%chebU(nn - 1)
      end do
    end if
  end subroutine chebTU

  subroutine simpBessel (this, r) 
    implicit none
    class(RadialFunction), intent(inout) :: this
    real(kind_double), intent(in) :: r
    integer :: kk, nn  
    real(kind_double) :: fact_kk, tmp_r, tmp_00, r_cut 
    
    r_cut = this%r_cut_out 
    tmp_00 = one_pi / r_cut
    tmp_r = r * tmp_00
    this%bg(0) = 1.d0
    this%bf(0) = 2.d0*one_pi/sqrt(2.5d0*r_cut**3) *(sincr(tmp_r) +sincr(2.d0*tmp_r))
    this%d_bf(0) = 2.d0*one_pi/sqrt(2.5d0*r_cut**3)*tmp_00 *(d_sincr(tmp_r) +2.d0*d_sincr(2.d0*tmp_r))

    do kk = 1, this%n_bessel
      fact_kk = dble((-1)**(kk+1) * (kk+1)*(kk+2))*one_pi/sqrt(dble((kk+1)**2 + (kk+2)**2)*0.5d0*r_cut**3)
      this%bf(kk) = fact_kk * (sincr(tmp_r*(kk+1) ) + sincr(tmp_r * (kk+2)))
      this%d_bf(kk) = fact_kk * tmp_00* ( (kk+1)*d_sincr(tmp_r*(kk+1)) + (kk+2)*sincr(tmp_r * (kk+2)) )

      this%bg(kk) = (this%bf(kk) + this%bg(kk-1)*sqrt(this%e_bess(kk)/this%d_bess(kk-1)))/sqrt(this%d_bess(kk))                         
      this%d_bg(kk) = (this%d_bf(kk) + this%d_bg(kk-1)*sqrt(this%e_bess(kk)/this%d_bess(kk-1)))/sqrt(this%d_bess(kk))                         
    end do

    do nn = 0, this%nradial
      ! nn = 0 is constant 1 so skipped in radial 
      this%only_radial(nn) = this%bg(nn)
      this%d_only_radial(nn) = this%d_bg(nn)
    end do 

  end subroutine simpBessel 

  subroutine radial_powTcheb (this, rr)
    implicit none
    class(RadialFunction), intent(inout) :: this
    !  for ACE the rbf goes from 0 to n_cheb-1
    real(kind_double), intent(in) :: rr
    real (kind_double) :: dtmp, tmp, xx, fderiv, r_cut 
    integer :: nn

    r_cut = this%r_cut_out 
    ! this%lambda should be > 1. 
    tmp = (1.d0 - rr / r_cut)
    dtmp = tmp ** (this%lambda-1)
    xx = 1.d0 - 2.d0 * dtmp*tmp


    fderiv = 2.d0*this%lambda / r_cut * dtmp

    call chebTU(this, xx)

    this%only_radial (0) = 1.d0 
    this%d_only_radial (0) = 0.d0 

    do nn = 1, this%n_cheb
      this%only_radial (nn) = 0.5d0 - 0.5d0 * this%chebT(nn)
      this%d_only_radial (nn) = -0.5d0 * fderiv * this%d_chebT(nn)
    end do

  return
end subroutine radial_powTcheb


subroutine radial_expTcheb (this, rr)
    implicit none
    class(RadialFunction), intent(inout) :: this
    !  for ACE the rbf goes from 0 to n_cheb-1
    real(kind_double), intent(in) :: rr
    real(kind_double) :: lambt, cost, d_cost, xx, d_xx, expt, tmpr, r_cut  
    integer :: nn    

    r_cut = this%r_cut_out 
    tmpr = rr / r_cut 

    cost = 1.d0 + cos(one_pi * tmpr )
    d_cost = -sin(one_pi * tmpr) * one_pi / r_cut 

    lambt = exp(this%lambda) - 1.d0 
    expt = exp(-this%lambda*(tmpr-1.d0))
    xx = 1.d0 - 2.d0*(expt  - 1.d0 ) / lambt 
    d_xx = 2.d0 * this%lambda / r_cut * expt / lambt 

    call chebTU(this, xx)

    this%only_radial(0) = 1.d0
    this%d_only_radial(0) = 0.d0 

    this%only_radial(1) = cost*0.5d0
    this%d_only_radial(1) = d_cost *0.5d0 

    do nn = 2, this%n_cheb
      this%only_radial(nn) = 0.25d0 * (1.d0 - this%chebT(nn-1)) * cost 
      this%d_only_radial(nn) = 0.25d0 * (1.d0 - this%chebT(nn-1)) *d_cost - 0.25d0 * d_xx*this%d_chebT(nn-1) * cost
    end do 

end subroutine radial_expTcheb 

function sincr(xx) result(tmp)
  implicit none
  real(kind_double), intent(in)  :: xx
  real(kind_double) ::  tmp

  if (xx == 0.d0 ) then
    tmp = 1.d0
  else
    tmp=sin(xx) / xx
  end if
end function sincr

function d_sincr(xx) result(tmp)
  implicit none
  real(kind_double), intent(in)  :: xx
  real(kind_double) ::  tmp

  if (xx == 0.d0 ) then
    tmp = 0.d0
  else
    tmp = (cos(xx) * xx - sin(xx)) / (xx**2)
  end if

end function d_sincr

function fcut(rcut, deltacut, rr) result(tmp)
  implicit none 
  real(kind_double), intent(in)  :: rcut, deltacut 
  real(kind_double) :: rr
  ! tmp(1) is the function, tmp(2) is the derivative. 
  real(kind_double), dimension(2) :: tmp

  real(kind_double) ::  Dtmp, yy, d_yy, xx, d_xx  
  integer :: val 

  Dtmp = rcut - deltacut 
 
  
  if (rr <= Dtmp ) then 
    tmp(1) = 1.d0 
    tmp(2) = 0.d0 
    val = 1 
  else if ( ( rr > Dtmp ).and.(rr < rcut) ) then 
    xx = 1.d0 - 2.d0 * ( 1.d0 + (rr - rcut)/deltacut )
    d_xx = -2.d0 / deltacut

    yy = 1.875d0 * xx - 1.25d0 * xx**3 + 0.375d0 * xx**5
    d_yy = (1.875d0 - 3.75d0 * xx**2  + 1.875d0 * xx**4 ) * d_xx 

    tmp(1) = 0.5d0 * (1.d0 + yy )   
    tmp(2) = 0.5d0 * d_yy 
    val = 2 

  else 
    tmp(1) = 0.d0 
    tmp(2) = 0.d0 
    val = 3 

  end if   

end function fcut 

subroutine evaluate_radial (this, rr)
  implicit none 
  class(RadialFunction), intent(inout) :: this
  real(kind_double) :: rr
  real(kind_double), dimension(2) :: cutout, cutin  

  cutout(:) = fcut(this%r_cut_out, this%delta_cut_out, rr) 
  cutin(:) = fcut(this%r_cut_in, this%delta_cut_in, rr) 

  if (this%radial_type == 'simpBessel' ) then 
    call simpBessel(this, rr)
    this%radial(:) = this%bg(0:this%n_bessel)  * cutout(1) * (1.d0 - cutin(1))
    this%d_radial(:) = this%d_bg(0:this%n_bessel)  * cutout(1) * (1.d0 - cutin(1)) + this%bg(:) * cutout(2) * (1.d0 - cutin(1)) - this%bg(:) * cutout(1) * cutin(2)

  else if (this%radial_type == 'powTcheb') then
    
      call radial_powTcheb(this, rr)
      this%radial(:) = this%only_radial(:)  * cutout(1) * (1.d0 - cutin(1))
      this%d_radial(:) = this%d_only_radial(:)  * cutout(1) * (1.d0 - cutin(1)) + this%only_radial(:) * cutout(2) * (1.d0 - cutin(1)) - this%only_radial(:) * cutout(1) * cutin(2) 

  else if (this%radial_type == 'expTcheb') then
      call radial_expTcheb(this, rr)

      this%radial(:) = this%only_radial(:)  * cutout(1) * (1.d0 - cutin(1))   
      this%d_radial(:) = this%d_only_radial(:)  * cutout(1) * (1.d0 - cutin(1)) + this%only_radial(:) * cutout(2) * (1.d0 - cutin(1)) - this%only_radial(:) * cutout(1) * cutin(2)

  else if (this%radial_type == 'powWave') then
      call log_critical('powWave not implemented yet.')
      
  else
      call log_critical('Unknown radial function type in ACE!')
      stop
  end if
    
end subroutine evaluate_radial  


end module RadialFunctions




module mod_radial_functions
  use module_kind_variables, only: kind_double
  implicit none
  private
  public :: paftouny_so3, sc_stefanodg

  contains


  subroutine paftouny_so3(desc_forces_local, n_cheb, xx)

    use module_so3, only: chebT, chebU, d_chebT

    logical, intent(in)  :: desc_forces_local
    real(kind_double), intent(in)    :: xx
    integer, intent(in)  :: n_cheb
    integer  :: nn

    chebT(0) = 1.d0
    chebT(1) = xx
    chebU(0) = 1.d0
    chebU(1) = 2.d0*xx

    d_chebT(0) = 0.d0
    d_chebT(1) = 1.d0
    if (n_cheb .ge. 2) then
      do nn = 2, n_cheb
        chebT(nn) = xx*chebT(nn - 1) - (1.d0 - xx**2)*chebU(nn - 2)
        chebU(nn) = xx*chebU(nn - 1) + chebT(nn)
      end do

      if (desc_forces_local) then
        do nn = 2, n_cheb
          d_chebT(nn) = dble(nn)*chebU(nn - 1)
        end do
      end if
    end if
  end subroutine paftouny_so3



  subroutine sc_stefanodg(desc_forces_local, n_cheb, rr)

    use ml_in_ndm_module, only: one_pi
    use module_neigh_local, only: r_cut 
    use module_so3, only: chebT, d_chebT, scgg, d_scgg

    integer, intent(in)  :: n_cheb
    logical, intent(in)  :: desc_forces_local
    real(kind_double), intent(in)    :: rr
    real(kind_double)    :: xx
    integer  :: nn

    xx = one_pi*(rr/r_cut)
    scgg(0) = 0.d0
    scgg(1) = 0.5d0*(1.d0 + dcos(xx))

    d_scgg(0) = 0.d0
    d_scgg(1) = -0.5d0*dsin(xx)*one_pi/r_cut

    if (n_cheb > 2) then
      do nn = 2, n_cheb
        scgg(nn) = 0.5d0*(1.d0 - chebT(nn - 1))*scgg(1)
        if (desc_forces_local) then
          ! THIS SHOULD BE x 2.D0/rcut
          d_scgg(nn) = -1.d0*scgg(1)*d_chebT(nn - 1)/r_cut + 0.5d0*(1.d0 - chebT(nn - 1))*d_scgg(1)
        end if
      end do
    end if
  end subroutine sc_stefanodg

end module mod_radial_functions
