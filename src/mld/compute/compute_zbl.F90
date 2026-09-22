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

#include "../../MLD_MACROS.INC"

module  module_potential_zbl 
  implicit none 
  contains 
  subroutine potential_zbl(desc_forces_local, Z_ja, Z_ia, r_ji, r_energy_zbl, r_force_zbl, l_r1_zbl, l_r2_zbl)
    use module_kind_variables, only: kind_double
    use ml_in_ndm_module, only: one_pi
    use module_zbl, only: epsilon0, p1_phi, p2_phi, p3_phi, e1_phi, e2_phi, e3_phi, & 
                          zz_exp, pz_phi, &
                          r1_zbl, r2_zbl 
#ifdef MLD_NDM
    use notperiod_mod
#endif 
    implicit none 
    logical, intent(in)  :: desc_forces_local
    real(kind_double), intent(in) :: Z_ia, Z_ja
    real(kind_double), intent(in) :: r_ji 
    real(kind_double), intent(out) :: r_energy_zbl 
    real(kind_double), intent(out) :: r_force_zbl   
    real(kind_double), intent(in), optional :: l_r1_zbl, l_r2_zbl
    real(kind_double) :: cst, tmp_a, tmp, tmp_phi, tmp_f, tmp_X, &
                         d_tmp, d_tmp_phi, d_tmp_f, force_zbl 
    real(kind_double) :: my_r1_zbl, my_r2_zbl
  
    ! use pair-specific values if provided, otherwise fall back to global
    my_r1_zbl = r1_zbl
    my_r2_zbl = r2_zbl
    if (present(l_r1_zbl)) my_r1_zbl = l_r1_zbl
    if (present(l_r2_zbl)) my_r2_zbl = l_r2_zbl
  
    cst = (Z_ia * Z_ja)/(4.d0 * one_pi * epsilon0)
    tmp_a = (Z_ia**zz_exp + Z_ja**zz_exp) / pz_phi !1/a
    tmp = r_ji * tmp_a !r_ij/a
    tmp_phi = p1_phi * dexp(e1_phi * tmp) + & 
              p2_phi * dexp(e2_phi * tmp) + &
              p3_phi * dexp(e3_phi * tmp)
  
  
    if (r_ji <= my_r1_zbl) then
      tmp_f = 1.d0
    else if (r_ji >= my_r2_zbl) then 
      tmp_f = 0.d0 
    else   
      tmp_X = (r_ji - my_r1_zbl)/(my_r2_zbl - my_r1_zbl)
      tmp_f = 1.d0 - tmp_X**3*(6.d0*tmp_X**2 - 15.d0*tmp_X + 10.d0)
    end if
  
    r_energy_zbl =  cst * tmp_phi * tmp_f / r_ji
  
    if (desc_forces_local) then
      d_tmp=tmp_a 
      d_tmp_phi = (p1_phi * e1_phi * dexp(e1_phi * tmp) + & 
                   p2_phi * e2_phi * dexp(e2_phi * tmp) + &
                   p3_phi * e3_phi * dexp(e3_phi * tmp))*d_tmp 
      if (r_ji <= my_r1_zbl) then
        d_tmp_f = 0.d0
      else if (r_ji >= my_r2_zbl ) then 
        d_tmp_f = 0.d0 
      else     
        d_tmp_f = (- 30.d0 * tmp_X**4 + 60.d0 * tmp_X**3 - 30.d0 * tmp_X**2)/(my_r2_zbl - my_r1_zbl)
      end if
      force_zbl = cst * ( (tmp_phi*d_tmp_f + d_tmp_phi*tmp_f)/r_ji - tmp_phi * tmp_f / r_ji**2)  !dV/dr_ij
      !flocal_zbl(1:3)  = force_zbl * dxp_ji(1:3)/r_ji/2.d0
      r_force_zbl   = force_zbl
    end if 
  
  end subroutine potential_zbl    
  
  subroutine potential_zbl_second(Z_ja, Z_ia, r_ji, r_energy_zbl, r_force_zbl, r_second_zbl, l_r1_zbl, l_r2_zbl)
    use module_kind_variables, only: kind_double
    use ml_in_ndm_module, only: one_pi
    use module_zbl, only: epsilon0, p1_phi, p2_phi, p3_phi, e1_phi, e2_phi, e3_phi, & 
                          zz_exp, pz_phi, &
                          r1_zbl, r2_zbl 
    implicit none 
    real(kind_double), intent(in) :: Z_ia, Z_ja
    real(kind_double), intent(in) :: r_ji 
    real(kind_double), intent(out) :: r_energy_zbl 
    real(kind_double), intent(out) :: r_force_zbl   
    real(kind_double), intent(out) :: r_second_zbl   
    real(kind_double), intent(in), optional :: l_r1_zbl, l_r2_zbl
    real(kind_double) :: cst, inv_a, inv_rij, xx,  chi, &
                        fcut_zbl, d_fcut_zbl, dd_fcut_zbl, &
                        phi, d_phi, dd_phi
    real(kind_double) :: my_r1_zbl, my_r2_zbl
  
    ! use pair-specific values if provided, otherwise fall back to global
    my_r1_zbl = r1_zbl
    my_r2_zbl = r2_zbl
    if (present(l_r1_zbl)) my_r1_zbl = l_r1_zbl
    if (present(l_r2_zbl)) my_r2_zbl = l_r2_zbl
  
    cst = (Z_ia * Z_ja)/(4.d0 * one_pi * epsilon0)
    inv_a = (Z_ia**zz_exp + Z_ja**zz_exp) / pz_phi !1/a
    xx = r_ji * inv_a !r_ij/a
  
    phi    =  p1_phi * dexp(e1_phi * xx) + &
              p2_phi * dexp(e2_phi * xx) + &
              p3_phi * dexp(e3_phi * xx)
    d_phi  = (p1_phi * e1_phi * dexp(e1_phi * xx) + & 
              p2_phi * e2_phi * dexp(e2_phi * xx) + &
              p3_phi * e3_phi * dexp(e3_phi * xx))*inv_a
    dd_phi = (p1_phi * e1_phi**2 * dexp(e1_phi * xx) + & 
              p2_phi * e2_phi**2 * dexp(e2_phi * xx) + &
              p3_phi * e3_phi**2 * dexp(e3_phi * xx))*inv_a**2
  
  
    if (r_ji <= my_r1_zbl) then
      fcut_zbl = 1.d0
      d_fcut_zbl = 0.d0
      dd_fcut_zbl = 0.d0 
    else if (r_ji >= my_r2_zbl) then 
      fcut_zbl = 0.d0 
      d_fcut_zbl = 0.d0 
      dd_fcut_zbl = 0.d0 
    else   
      chi = (r_ji - my_r1_zbl)/(my_r2_zbl - my_r1_zbl)
      fcut_zbl = 1.d0 - chi**3*(6.d0*chi**2 - 15.d0*chi + 10.d0)
      d_fcut_zbl = (- 30.d0 * chi**4 + 60.d0 * chi**3 - 30.d0 * chi**2)/(my_r2_zbl - my_r1_zbl)
      dd_fcut_zbl = (- 120.d0 * chi**3 + 180.d0 * chi**2 - 60.d0 * chi)/(my_r2_zbl - my_r1_zbl)**2
    end if
    inv_rij = 1.d0 /  r_ji
    r_energy_zbl =  cst * phi * fcut_zbl * inv_rij
    r_force_zbl = cst * inv_rij * ( phi*d_fcut_zbl + d_phi*fcut_zbl  - phi * fcut_zbl * inv_rij )  !dV/dr_ij
    r_second_zbl =    cst * inv_rij * (phi * dd_fcut_zbl + 2.d0 * d_phi * d_fcut_zbl + dd_phi * fcut_zbl )  &
                    - cst * inv_rij**2 * ( d_phi * fcut_zbl + phi * d_fcut_zbl)  &
                    + cst * 2.d0 * inv_rij**3 * phi * fcut_zbl   
     
  
  end subroutine potential_zbl_second    
  
  
  subroutine full_potential_zbl (desc_forces_local, spair_potential, r_ji, r_energy_zbl, r_force_zbl)
    use module_kind_variables, only: kind_double
    use module_zbl, only: params_k2b_to_zbl, r1_zbl, rr_k2b, type_rac_zbl, rac_zbl_exp, rac_zbl_poly, &
                          r1_zbl_pair, rr_k2b_pair
    use derived_types, only: typ_species_half
  
    implicit none 
    integer, intent(in)  :: spair_potential 
    logical, intent(in)  :: desc_forces_local 
    real(kind_double), intent(in) :: r_ji 
    real(kind_double), intent(out) :: r_energy_zbl, r_force_zbl 
    real(kind_double) :: Z_ia, Z_ja
    real(kind_double) :: my_r1_zbl, my_rr_k2b
    integer :: itype1, itype2
  
    Z_ja = dble(typ_species_half(spair_potential)%Z1)
    Z_ia = dble(typ_species_half(spair_potential)%Z2)
    itype1 = typ_species_half(spair_potential)%type1
    itype2 = typ_species_half(spair_potential)%type2

    ! use pair-specific values if available
    if (allocated(r1_zbl_pair)) then
      my_r1_zbl = r1_zbl_pair(itype1, itype2)
      my_rr_k2b = rr_k2b_pair(itype1, itype2)
    else
      my_r1_zbl = r1_zbl
      my_rr_k2b = rr_k2b
    end if
  
    if (r_ji <= my_r1_zbl) then
      call potential_zbl(desc_forces_local, Z_ja, Z_ia, r_ji, r_energy_zbl, r_force_zbl, &
                         l_r1_zbl=my_r1_zbl) 
    else if ((r_ji <= my_rr_k2b).and.(r_ji > my_r1_zbl)) then 
      select case( type_rac_zbl )
        
        case ( rac_zbl_exp )  
          r_energy_zbl = dexp( params_k2b_to_zbl(1,spair_potential)           +   &
                           params_k2b_to_zbl(2,spair_potential)*r_ji      +   &
                           params_k2b_to_zbl(3,spair_potential)*r_ji**2   +   &
                           params_k2b_to_zbl(4,spair_potential)*r_ji**3)
          if (desc_forces_local) then 
             r_force_zbl = (     params_k2b_to_zbl(2,spair_potential)          +  &
                       2.d0*params_k2b_to_zbl(3,spair_potential)*r_ji     +  &
                       3.d0*params_k2b_to_zbl(4,spair_potential)*r_ji**2) *  &
                       r_energy_zbl                 
          end if 
        case ( rac_zbl_poly ) 
          r_energy_zbl = params_k2b_to_zbl(1,spair_potential)           +   &
                     params_k2b_to_zbl(2,spair_potential)*r_ji      +   &
                     params_k2b_to_zbl(3,spair_potential)*r_ji**2   +   &
                     params_k2b_to_zbl(4,spair_potential)*r_ji**3   +   &
                     params_k2b_to_zbl(5,spair_potential)*r_ji**4   +   &
                     params_k2b_to_zbl(6,spair_potential)*r_ji**5
          if (desc_forces_local) then 
            r_force_zbl =     params_k2b_to_zbl(2,spair_potential)           +   &
                     2.d0*params_k2b_to_zbl(3,spair_potential)*r_ji      +   &
                     3.d0*params_k2b_to_zbl(4,spair_potential)*r_ji**2   +   &
                     4.d0*params_k2b_to_zbl(5,spair_potential)*r_ji**3   +   &
                     5.d0*params_k2b_to_zbl(6,spair_potential)*r_ji**4
          end if                  
      end select    
    else
      r_energy_zbl = 0.d0 
      r_force_zbl =  0.d0                        
    end if   
    end subroutine full_potential_zbl 

    subroutine zbl_end(Z_ja, Z_ia, r1, f1, f1d, f1dd)
      !
      ! provide the ZBL function f1 and derivative f1d at the distance r1 
      !
      use module_kind_variables, only: kind_double
      use mld_logger
      use module_zbl, only: type_rac_zbl, rac_zbl_exp, rac_zbl_poly
      implicit none 
      real(kind_double), intent(in)  :: r1, Z_ja, Z_ia 
      real(kind_double), intent(out) :: f1, f1d, f1dd 
  
      logical :: desc_forces_local 
      
      _NAMECURRENT_("zbl_end")
  
      !ckeck consistentcy itype1, itype2, and Z_ia, Z_ja !!!  
  
      desc_forces_local = .true. 
  
      !check id f1 and f1d is real the right one !!!!
      f1dd = 0.d0 
      select case (type_rac_zbl)

        case (rac_zbl_exp)
          call potential_zbl(desc_forces_local, Z_ja, Z_ia, r1,  f1, f1d)
        case (rac_zbl_poly)
          call potential_zbl_second (Z_ja, Z_ia, r1,  f1, f1d, f1dd )
      end select     

      _MLD_END_
  
    end subroutine zbl_end  
  


  end  module module_potential_zbl 
  

module module_compute_zbl 

contains 

subroutine compute_zbl(i_start_at, i_final_at, d_n_neigh, d_kind_neigh,  iconf)

#if(PARA)
  !use mpi
  use mld_mpi 
#endif
  USE module_kind_variables, ONLY: kind_double
#ifdef MLD_NDM
  use gen_com_m, ONLY: A2cm, lperiod
  use gen_com_m_ml, ONLY: imm, bg, at
  use tab_imm_m_ml, ONLY: xp
#else
  use ondm_gen_com_m, ONLY: imm, lperiod
  use ondm_tab_imm_m, ONLY: iwmax2, xp
#endif
  use derived_types, only: config_real
  use ml_in_ndm_module, ONLY: rangml, imm_neigh, desc_forces
  use module_neigh_local, only: r_central, i_central, i_type, i_type_db, tmp_dxp, tmp_xp, &
                                max_neigh_local, iw2, build_local_neighbours_ja
  use time_check_general, only: MY_MPI_WTIME
  use module_zbl, only: r2_zbl, r1_zbl_pair, r2_zbl_pair
  use module_potential_zbl, only: potential_zbl                      
  use mld_logger
#ifdef MLD_NDM
  use notperiod_mod
#else
  use ondm_transform_coord, only: ondm_notperiod
#endif
  implicit none

  integer, intent(in)  :: i_start_at, i_final_at
  integer, dimension(imm), intent(out)   :: d_n_neigh
  integer, dimension(imm, imm_neigh), intent(out)    :: d_kind_neigh
  integer, optional    :: iconf
  logical  :: small
  real(kind_double), dimension(:, :), allocatable   :: xpnp
  real(kind_double), dimension(3) :: dxp_ji
  integer  :: ia, ja, ia_n
  real(kind_double)   :: r_ji
  logical  :: desc_forces_local
  real(kind_double)    :: Z_ia, Z_ja
  real(kind_double)    :: energy_zbl, force_zbl, r_energy_zbl
  real(kind_double), dimension(6)     :: stress_zbl
  integer :: itype_ja, itype_ia
  real(kind_double) :: my_r1_zbl, my_r2_zbl
  _NAMECURRENT_("compute_zbl") 

  _MLD_BEGIN_

  !if (not(allocated(config_real(iconf)%elocal_zbl))) allocate(config_real(iconf)%elocal_zbl(config_real(iconf)%nat)) 
  config_real(iconf)%elocal_zbl(:) = 0.d0

  if ((i_start_at == 0) .and. (i_final_at == 0)) then
    d_n_neigh(:) = 0
    d_kind_neigh(:, :) = 0
    return
  end if

  small = .false.
  if (present(iconf)) then
    small = config_real(iconf)%small
  end if


  desc_forces_local = desc_forces .and. (config_real(iconf)%has_force .or. config_real(iconf)%has_stress)
  ALLOCATE (xpnp(3, imm))
  if (lperiod) then
    xpnp(:, :) = xp(:, :)
  else
#ifdef MLD_NDM
         call notperiod(imm,xp, xpnp,at,bg,.false.)
#else
    call ondm_notperiod(xp, xpnp)
#endif
  end if

  d_n_neigh(:) = 0
  d_kind_neigh(:, :) = 0
  config_real(iconf)%elocal_zbl(:) = 0.d0
  if (desc_forces_local) then 
    !if (allocated(config_real(iconf)%flocal_zbl))  deallocate(config_real(iconf)%flocal_zbl)
    config_real(iconf)%flocal_zbl(:, :, :) = 0.d0
  end if 

  if (config_real(iconf)%has_stress) then
    if (allocated(config_real(iconf)%slocal_zbl))  deallocate(config_real(iconf)%slocal_zbl)
    allocate(config_real(iconf)%slocal_zbl(config_real(iconf)%nat,6))
    config_real(iconf)%slocal_zbl(:, :) = 0.d0
  end if
#ifdef MLD_NDM
  iw2=0
!!$  if (i_start_at == 1) iw2 = 0
!!$  if (i_start_at > 1) iw2 = iwmax2(i_start_at - 1)
#else
  if (i_start_at == 1) iw2 = 0
  if (i_start_at > 1) iw2 = iwmax2(i_start_at - 1)
#endif

  do ja = i_start_at, i_final_at
    ! begin small box or not 1/

    call build_local_neighbours_ja( iconf, ja, imm, xpnp, r2_zbl, d_n_neigh, d_kind_neigh, &
                                    r_central, i_type, i_type_db, i_central, tmp_dxp, tmp_xp, &
                                    max_neigh_local, iw2)

    energy_zbl = 0.d0 
    stress_zbl(:) = 0.d0
    Z_ja = config_real(iconf)%Z_per_type(config_real(iconf)%itype(ja))
    itype_ja = config_real(iconf)%itype(ja)

    do ia_n = 1, max_neigh_local
      ! begin small box or not 2/
      r_ji = r_central(ia_n)
      !write(*,*) rangml, 'r_ji', r_ji, ia_n, max_neigh_local

      ia = i_central(ia_n)      
      dxp_ji = tmp_dxp(1:3, ia_n)
      Z_ia = config_real(iconf)%Z_per_type(config_real(iconf)%itype(ia))
      itype_ia = config_real(iconf)%itype(ia)

      ! lookup pair-specific ZBL boundaries
      if (allocated(r1_zbl_pair)) then
        my_r1_zbl = r1_zbl_pair(itype_ja, itype_ia)
        my_r2_zbl = r2_zbl_pair(itype_ja, itype_ia)
      else
        my_r1_zbl = 0.d0
        my_r2_zbl = r2_zbl
      end if

      call potential_zbl(desc_forces_local, Z_ja, Z_ia, r_ji, r_energy_zbl, force_zbl, &
                         l_r1_zbl=my_r1_zbl, l_r2_zbl=my_r2_zbl)
      energy_zbl = energy_zbl + r_energy_zbl 
      if (desc_forces_local) then 
        config_real(iconf)%flocal_zbl(ja, ia_n, 1:3)  = force_zbl * dxp_ji(1:3)/r_ji/2.d0
      end if 
      

      if (config_real(iconf)%has_stress) then
        stress_zbl(1) = stress_zbl(1) + config_real(iconf)%flocal_zbl(ja, ia_n, 1)*dxp_ji(1)/r_ji
        stress_zbl(2) = stress_zbl(2) + config_real(iconf)%flocal_zbl(ja, ia_n, 2)*dxp_ji(2)/r_ji
        stress_zbl(3) = stress_zbl(3) + config_real(iconf)%flocal_zbl(ja, ia_n, 3)*dxp_ji(3)/r_ji
        stress_zbl(4) = stress_zbl(4) + config_real(iconf)%flocal_zbl(ja, ia_n, 2)*dxp_ji(3)/r_ji
        stress_zbl(5) = stress_zbl(5) + config_real(iconf)%flocal_zbl(ja, ia_n, 1)*dxp_ji(3)/r_ji
        stress_zbl(6) = stress_zbl(6) + config_real(iconf)%flocal_zbl(ja, ia_n, 1)*dxp_ji(2)/r_ji
      end if

    end do       ! ia end of neighbours iterations

    ! energy final .......
    config_real(iconf)%elocal_zbl(ja) = energy_zbl/2.d0

    ! forces final .......
    if (desc_forces_local) then
      if (.not.((i_start_at == 0) .and. (i_final_at == 0))) then 
      config_real(iconf)%flocal_zbl(ja, 0, 1:3)  = - sum(config_real(iconf)%flocal_zbl(ja, 1:max_neigh_local,1:3), dim=1)
      end if 
    end if

    ! stress final .......
    if (config_real(iconf)%has_stress) then
      config_real(iconf)%slocal_zbl(ja, 1) = stress_zbl(1)
      config_real(iconf)%slocal_zbl(ja, 2) = stress_zbl(2)
      config_real(iconf)%slocal_zbl(ja, 3) = stress_zbl(3)
      config_real(iconf)%slocal_zbl(ja, 4) = stress_zbl(4)
      config_real(iconf)%slocal_zbl(ja, 5) = stress_zbl(5)
      config_real(iconf)%slocal_zbl(ja, 6) = stress_zbl(6)
    end if
  end do

  _MLD_END_
end subroutine compute_zbl



subroutine energy_force_stress_zbl (iconf)
#if (PARA)
  use mpi 
  use mld_mpi 
#endif
#ifdef MLD_NDM  
  use gen_com_m_ml, only : imm
#else
 use ondm_gen_com_m, only : imm
#endif
  use module_kind_variables, only: kind_double
  use derived_types, only: config_real
  use ml_in_ndm_module, only: rangml, i_start_at, i_final_at, desc_forces, &
                              imm_neigh
  use mld_logger
  use mld_subworld, only: subworld
  use my_mpi_subroutines, only: my_barrier_subworld
  implicit none 
  integer, intent(in) :: iconf 
  integer :: dim_reduce, ik, inn1, ia
  logical :: desc_forces_local 
  _NAMECURRENT_("energy_force_stress_zbl") 

  _MLD_BEGIN_

#if (PARA)


  !call MPI_BARRIER(subworld, codeml)
  call my_barrier_subworld(codeml)
  call MPI_ALLREDUCE(MPI_IN_PLACE, config_real(iconf)%n_neigh_zbl, imm, MPI_INTEGER, MPI_SUM, subworld, codeml)
  call MPI_ALLREDUCE(MPI_IN_PLACE, config_real(iconf)%kind_neigh_zbl, imm*imm_neigh, MPI_INTEGER, MPI_SUM, subworld, codeml)
  !call MPI_BARRIER(subworld, codeml)
  call my_barrier_subworld(codeml)

#endif

  desc_forces_local = desc_forces .and. (config_real(iconf)%has_force .or. config_real(iconf)%has_stress)

  !zbl energy 
  dim_reduce = config_real(iconf)%nat
  call MPI_ALLREDUCE(MPI_IN_PLACE, config_real(iconf)%elocal_zbl, dim_reduce, MPI_DOUBLE_PRECISION, MPI_SUM, subworld, codeml)
  config_real(iconf)%ezbl = sum(config_real(iconf)%elocal_zbl(:))
  if (allocated(config_real(iconf)%elocal_zbl)) deallocate(config_real(iconf)%elocal_zbl) 


  if (desc_forces_local) then 
    !zbl forces
    if (allocated(config_real(iconf)%fzbl)) deallocate(config_real(iconf)%fzbl)
    allocate(config_real(iconf)%fzbl(3,config_real(iconf)%nat))
    config_real(iconf)%fzbl(:,:) = 0.d0  
    
    if (.not.((i_start_at == 0) .and. (i_final_at == 0))) then 
      do ik = i_start_at, i_final_at                   ! derivative index \partial / \partial ik,alpha
        config_real(iconf)%fzbl(1:3, ik) = config_real(iconf)%fzbl(1:3, ik) - config_real(iconf)%flocal_zbl(ik, 0, 1:3)
        ! icount1=0
        do inn1 = 1, config_real(iconf)%n_neigh_zbl(ik)      ! sum over all atoms
          ia = config_real(iconf)%kind_neigh_zbl(ik, inn1)
          config_real(iconf)%fzbl(1:3, ia) = config_real(iconf)%fzbl(1:3, ia) - config_real(iconf)%flocal_zbl(ik, inn1, 1:3)
        end do
      end do
    end if 

    dim_reduce = 3*config_real(iconf)%nat
    call MPI_ALLREDUCE(MPI_IN_PLACE, config_real(iconf)%fzbl, dim_reduce, MPI_DOUBLE_PRECISION, MPI_SUM, subworld, codeml)

    if (allocated(config_real(iconf)%flocal_zbl)) deallocate(config_real(iconf)%flocal_zbl)
  end if


  ! zbl stress 
  
  if (config_real(iconf)%has_stress) then 

    if (allocated(config_real(iconf)%szbl)) deallocate(config_real(iconf)%szbl)
    allocate(config_real(iconf)%szbl(6))
    dim_reduce = 6*size(config_real(iconf)%slocal_zbl,dim=1)
    call MPI_ALLREDUCE(MPI_IN_PLACE, config_real(iconf)%slocal_zbl, dim_reduce, MPI_DOUBLE_PRECISION, MPI_SUM, subworld, codeml)
  
    config_real(iconf)%szbl(:) = sum(config_real(iconf)%slocal_zbl(:,:), dim=1)

    if (allocated(config_real(iconf)%slocal_zbl)) deallocate(config_real(iconf)%slocal_zbl)
  end if
   
 
  _MLD_END_
end subroutine energy_force_stress_zbl

subroutine allocate_zbl(iconf, i_start_at, i_final_at,  desc_forces_local)
#ifdef MLD_NDM
  use ml_in_ndm_module, only : imm_neigh, rangml
  use gen_com_m_ml, only:imm
#else
  use ml_in_ndm_module, only : imm_neigh, imm, rangml
#endif
  use derived_types, only: config_real 
  use mld_logger
  implicit none 
  integer, intent(in) :: iconf 
  integer, intent(in) :: i_start_at, i_final_at
  logical, intent(in) :: desc_forces_local 

  _NAMECURRENT_("allocate_zbl")
  if (imm /= config_real(iconf)%nat) then
    write(*,*) rangml, 'error in allocate_zbl'
    stop
  end if
  !energy 
  if (allocated(config_real(iconf)%elocal_zbl)) deallocate (config_real(iconf)%elocal_zbl) 
  allocate (config_real(iconf)%elocal_zbl(imm))

  config_real(iconf)%elocal_zbl = 0.d0

  !force 
  if (desc_forces_local) then
    if (.not. ((i_start_at == 0) .and. (i_final_at == 0))) then
      if (allocated(config_real(iconf)%flocal_zbl)) deallocate (config_real(iconf)%flocal_zbl) 
      allocate (config_real(iconf)%flocal_zbl(i_start_at:i_final_at, 0:imm_neigh, 1:3))
    end if
  end if

  if (allocated(config_real(iconf)%flocal_zbl)) then 
   config_real(iconf)%flocal_zbl = 0.d0 
  end if

  !stress
  if (config_real(iconf)%has_stress) then
    if (allocated(config_real(iconf)%slocal_zbl))  deallocate(config_real(iconf)%slocal_zbl)
    allocate(config_real(iconf)%slocal_zbl(config_real(iconf)%nat,6))
  end if

  if (allocated(config_real(iconf)%slocal_zbl)) then 
   config_real(iconf)%slocal_zbl = 0.d0 
  end if 

  _MLD_END_
end subroutine allocate_zbl 

end module module_compute_zbl 


module module_nlinear_compute
  use temporary_data_cov, only: dim_xdesc
  implicit none 
  contains  

  subroutine init_nlinear
    use module_nlinear, only : alpha_nl, beta_nl, l_nl_order,nl_order, &
                                nparams_nlinear
    
    implicit none 

    nparams_nlinear = 0 
    if (l_nl_order(2)) then 
      if (allocated(alpha_nl)) deallocate (alpha_nl)
      allocate (alpha_nl(nl_order(2), dim_xdesc))
      nparams_nlinear = nparams_nlinear + nl_order(2)*dim_xdesc 
    end if 
     
    if (l_nl_order(3)) then 
      if (allocated(beta_nl)) deallocate (beta_nl)
      allocate (beta_nl(nl_order(3), dim_xdesc))
      nparams_nlinear = nparams_nlinear + nl_order(3)*dim_xdesc 
    end if 
  
  end subroutine init_nlinear

  subroutine evaluate_nlinear_2 (desc, d_desc, max_neigh_local, desc_forces, e_tmp, f_tmp, &
            lJweig, e_weig_tmp, f_weig_tmp)
  ! TODO 
  ! dimension descc ans ddesc are nor logical ... 
  use module_kind_variables, only: kind_double
  use module_nlinear, only : alpha_nl 
  implicit none 
  real(kind_double), dimension(:), intent(in) :: desc    
  real(kind_double), dimension(:,:,:), intent(in) :: d_desc 
  logical, intent(in) :: desc_forces, lJweig    
  integer, intent(in) :: max_neigh_local 
  real(kind_double), intent(out) :: e_tmp
  real(kind_double), dimension(:,:), intent(out) :: f_tmp
  real(kind_double), dimension(:,:), intent(out) :: e_weig_tmp
  real(kind_double), dimension(:,:,:,:), intent(out) :: f_weig_tmp
  real(kind_double) :: tmp1, tmp2, dtmp1, dtmp2  
  integer  :: inn, ix, id 
  real(kind_double), dimension(:,:), allocatable :: f1, f2 
  real(kind_double), dimension(:), allocatable :: ftmp1, ftmp2 

    allocate (f1(1:dim_xdesc,3))
    allocate (f2(1:dim_xdesc,3))
    allocate(ftmp1(1:dim_xdesc))
    allocate(ftmp2(1:dim_xdesc))

    e_tmp = 0.d0
    tmp1 = dot_product(alpha_nl(1, 1:dim_xdesc),desc(2:dim_xdesc+1) )
    tmp2 = dot_product(alpha_nl(2, 1:dim_xdesc),desc(2:dim_xdesc+1) )
    e_tmp =  tmp1*tmp2 
    if (lJweig) then
      do id = 1, dim_xdesc
        e_weig_tmp(id, 1) = tmp2 * desc(1 + id) 
        e_weig_tmp(id, 1) = tmp1 * desc(1 + id) 
      end do 
    end if 
    if (desc_forces) then 
      f1(:,:) = 0.d0  
      f2(:,:) = 0.d0  
      do inn = 1 , max_neigh_local
        do ix=1,3
           dtmp1= dot_product (alpha_nl(1, 1:dim_xdesc) , d_desc(1:dim_xdesc,inn,ix))
           dtmp2= dot_product (alpha_nl(2, 1:dim_xdesc) , d_desc(1:dim_xdesc,inn,ix))
           f_tmp(inn, ix) = tmp1 * dtmp2 + dtmp1 * tmp2
           if (lJweig) then
            ftmp1(1:dim_xdesc) = desc(2:dim_xdesc+1) * dtmp2 + d_desc(1:dim_xdesc,inn,ix) * tmp2  
            f_weig_tmp(1:dim_xdesc,inn,ix,1) = ftmp1(1:dim_xdesc)
            ftmp2(1:dim_xdesc)  = desc(2:dim_xdesc+1) * dtmp1 + d_desc(1:dim_xdesc,inn,ix) * tmp1  
            f_weig_tmp(1:dim_xdesc,inn,ix,2) = ftmp2(1:dim_xdesc)
            f1(1:dim_xdesc,ix) = f1(1:dim_xdesc,ix) - ftmp1(1:dim_xdesc)  
            f2(1:dim_xdesc,ix) = f2(1:dim_xdesc,ix) - ftmp2(1:dim_xdesc)  
           end if 
        end do 
      end do !inn 
      !TODO f_weig_tmp(1:dim_xdesc,0,1:3,2) = f2(1:dim_xdesc,1:3)
      !TODO f_weig_tmp(1:dim_xdesc,0,1:3,1) = f1(1:dim_xdesc,1:3)
    end if 

    deallocate(f1, f2, ftmp1, ftmp2)
  
  end subroutine evaluate_nlinear_2 


  subroutine evaluate_nlinear_3 (desc, d_desc, max_neigh_local, desc_forces, e_tmp, f_tmp)
    use module_kind_variables, only: kind_double
    use module_nlinear, only : beta_nl 
    implicit none 
    real(kind_double), dimension(:), intent(in) :: desc    
    real(kind_double), dimension(:,:,:), intent(in) :: d_desc 
    logical, intent(in) :: desc_forces    
    integer, intent(in) :: max_neigh_local 
    real(kind_double), intent(out) :: e_tmp
    real(kind_double), dimension(:,:), intent(out) :: f_tmp
    real(kind_double) :: tmp1, tmp2, tmp3, dtmp1, dtmp2, dtmp3  
    integer  :: inn, ix
  
      e_tmp = 0.d0
      tmp1 = dot_product(beta_nl(1, 1:dim_xdesc),desc(2:dim_xdesc+1) )
      tmp2 = dot_product(beta_nl(2, 1:dim_xdesc),desc(2:dim_xdesc+1) )
      tmp3 = dot_product(beta_nl(3, 1:dim_xdesc),desc(2:dim_xdesc+1) )
      e_tmp =  tmp1*tmp2*tmp3
      if (desc_forces) then 
        do inn = 1 , max_neigh_local
          do ix=1,3
             dtmp1= dot_product (beta_nl(1, 1:dim_xdesc) , d_desc(1:dim_xdesc,inn,ix))
             dtmp2= dot_product (beta_nl(2, 1:dim_xdesc) , d_desc(1:dim_xdesc,inn,ix))
             dtmp3= dot_product (beta_nl(3, 1:dim_xdesc) , d_desc(1:dim_xdesc,inn,ix))
             f_tmp(inn, ix) = dtmp1 * tmp2 * tmp3 + tmp1 * dtmp2 * tmp3  + tmp1 * tmp2 * dtmp3 
          end do 
        end do 
      end if 
    
    end subroutine evaluate_nlinear_3 
  


  end module module_nlinear_compute 
