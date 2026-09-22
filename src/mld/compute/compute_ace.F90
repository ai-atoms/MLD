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





module module_compute_ace
contains 
subroutine compute_ace(i_start_at, i_final_at, d_n_neigh, d_kind_neigh, iconf)
  use iso_fortran_env, only: dp => real64
#ifdef MLD_NDM  
  use gen_com_m, ONLY: A2cm, lperiod
  use gen_com_m_ml, ONLY: imm, bg, at
  use tab_imm_m_ml, ONLY: xp ! JPC je pense que xp ne sert à rien
#else
  use ondm_gen_com_m, ONLY: imm,  lperiod
  use ondm_tab_imm_m, ONLY: iwmax2, xp 
#endif
  use mld_logger
  use angular_functions, only: spherical_harm, grad_spherical_harm
  use derived_types, only: config_real, config_desc
  use ml_in_ndm_module, ONLY: imm_neigh, &
                              desc_forces, linvisible, weighted
  use module_ace_desc, only: l_ace_order, ace_numax, time_amat, time_amat_01, time_amat_02, time_amat_03, time_cgord, & 
                             time_amat_00, &
                             time_ylm, time_rad, time_unique, &
                             zetaace_order, delta_zetaace, pos_ace_chem, ace_kmax
  use module_neigh_local, only: r_central, i_central, i_type, i_type_db, tmp_dxp, tmp_xp, &
                                max_neigh_local, iw2, & 
                                build_local_neighbours_ja, reallocate_neigh_ja, build_neigh_ja_type
  use time_check_general, only: time, tot_time, debug_time, MY_MPI_WTIME
  !use module_so3, only: n_rbf_so3, clmn,  ini_rbf_so3
  use compute_pow_so3_mod, only: spherical_3d
  use module_neigh_local, only: r_cut 
  use module_AB_basis, only: compute_B_basis_atom, compute_ylm_for_neigh_ja, compute_rad_for_neigh_ja, &
                             compute_uniqueA_neigh_ja, compute_B_basis_atom_adjoint, compute_AB_body_one_neigh_ja
  use module_base_cnlm, only: base_cnlm
  
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
  
  real(dp), dimension(:), allocatable :: baseB
  real(dp), dimension(:,:,:), allocatable :: d_baseB


  logical  :: small
  real(dp), dimension(:, :), allocatable   :: xpnp
  real(dp) :: temp_dja, time_aa, time_bb 
  integer  :: ja, ix, ia_n, dim_baseB, icount, inu, ii , iz 

  !real(double), dimension(n_rbf_so3)     :: dphi_ji
  logical  :: desc_forces_local
  integer  :: ja_atom, type_db_ja, max_neigh, ja_pos_in_desc
  _NAMECURRENT_("compute_ace")
  _MLD_BEGIN_

  max_neigh = imm_neigh 
  if ((i_start_at == 0) .and. (i_final_at == 0)) then
    d_n_neigh(:) = 0
    d_kind_neigh(:, :) = 0
    config_desc(iconf)%energy(:, :) = 0.d0
    ! config_desc(iconf)%force(:,:,:,:)=0.d0
    return
  end if

  small = .false.
  if (present(iconf)) then
    small = config_real(iconf)%small
  end if

  desc_forces_local = desc_forces .and. (config_real(iconf)%has_force .or. config_real(iconf)%has_stress)

!JPC XPNP NE SERT A RIEN (inutile dans  build_local_neighbours_ja)
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
  config_desc(iconf)%energy(:, :) = 0.d0
  if (desc_forces_local) config_desc(iconf)%force(:, :, :, :) = 0.d0

  if (i_start_at == 1) iw2 = 0
#ifdef MLD_NDM
!JPC IW2 est l'indice du preimer voisin ca vaut toujours 1 pour small
!  if (i_start_at > 1) iw2 = iwmax2(i_start_at - 1)
  iw2=0  
!JPC 
#else
  if (i_start_at > 1) iw2 = iwmax2(i_start_at - 1)
#endif  

  if (debug_time) time(1) = MY_MPI_WTIME()
  do ja = i_start_at, i_final_at
    ! begin small box or not 1/
    ja_atom = ja
    type_db_ja = config_real(iconf)%itype_db(ja)
    !ACE_CHEM
    ja_pos_in_desc = (type_db_ja-1)*pos_ace_chem

    if (debug_time) time(3) = MY_MPI_WTIME()
    if (linvisible .and. weighted) then
      if (config_real(iconf)%invisible_per_type(config_real(iconf)%itype(ja))) cycle
    end if


    ia_n = 0
    call build_local_neighbours_ja( iconf, ja, imm, xpnp, r_cut, d_n_neigh, d_kind_neigh, &
                                    r_central, i_type, i_type_db, i_central, tmp_dxp, tmp_xp, &
                                    max_neigh_local, iw2)

    if (max_neigh_local == 0) cycle
    if (max_neigh_local > max_neigh) then
      call log_critical("ML: Fatal error in "//NAMECURRENT//" : number of neighbours exceeded.")
      call log_critical("max_neigh_local="// vtoa(max_neigh_local) //" > max_neigh=imm_neigh="// vtoa(max_neigh))
      call log_critical("config file: "//trim(config_real(iconf)%filename))
      call log_critical("Solutions: 1) decrease r_cut  2) increase imm_neigh in ml_in_ndm_module.F90") 
      stop 'compute_ace: max_neigh_local > imm_neigh. Increase imm_neigh in ml_in_ndm_module.F90'
    end if
    call reallocate_neigh_ja (max_neigh_local, ja_atom, i_central, i_type, i_type_db, r_central, tmp_dxp)
    call build_neigh_ja_type (max_neigh_local, ja_atom, i_central, i_type, i_type_db, r_central, tmp_dxp)
    !debug write(565, *) max_neigh_local, trim(config_real(iconf)%filename)
    
    if (debug_time) then
      time(4) = MY_MPI_WTIME()
      tot_time(2) = tot_time(2) + time(4) - time(3)
    end if
    
    icount = 0 
    if (l_ace_order(1)) then
          inu = 1 
          ! rsp/drsp are no longer used: radials are evaluated directly 
          ! in compute_AB_body_one_neigh_ja (via radial_spline or chem_compressor)
          !was! call  compute_rad_for_neigh_ja(inu)
          call  compute_AB_body_one_neigh_ja(dim_baseB, baseB, d_baseB, desc_forces_local)
          !ACE_CHEM
          config_desc(iconf)%energy(icount+1+ja_pos_in_desc : icount+dim_baseB+ja_pos_in_desc, ja) = baseB(1:dim_baseB)
          !config_desc(iconf)%energy(icount+1:icount+dim_baseB, ja) = baseB(1:dim_baseB)
          if (desc_forces_local) then 
            do ix = 1,3
              !ACE_CHEM
              config_desc(iconf)%force(icount+1+ja_pos_in_desc : icount+dim_baseB+ja_pos_in_desc, ja, 1:max_neigh_local, ix) = d_baseB(1:dim_baseB, ix, 1:max_neigh_local) 
              !config_desc(iconf)%force(icount+1:icount+dim_baseB, ja, 1:max_neigh_local, ix) = d_baseB(1:dim_baseB, ix, 1:max_neigh_local)
            end do        
          end if 
          icount = icount + dim_baseB

          
          do iz = 2, zetaace_order
            do ii = 1, dim_baseB
              icount = icount  + 1
              !$! temp_dja = ( delta_zetaace(iz,1) * baseB(ii)  )**(iz-1)
              !$! config_desc(iconf)%energy(icount, ja) =  temp_dja *  delta_zetaace(iz,1) * baseB(ii)
              temp_dja = delta_zetaace(iz,1) * ( baseB(ii) )**(iz-1)
              !ACE_CHEM
              config_desc(iconf)%energy(icount+ja_pos_in_desc, ja) =  temp_dja * baseB(ii)
              !config_desc(iconf)%energy(icount, ja) =  temp_dja * baseB(ii)
            
              if (desc_forces_local) then 
                do ix = 1, 3
                  !ACE_CHEM
                  config_desc(iconf)%force(icount+ja_pos_in_desc, ja, 1:max_neigh_local, ix) =  dble(iz) * d_baseB(ii, ix, 1:max_neigh_local) * temp_dja 
                  !config_desc(iconf)%force(icount, ja, 1:max_neigh_local, ix) =  dble(iz) * d_baseB(ii, ix, 1:max_neigh_local) * temp_dja 
                end do
              end if         
          end do 

        end do 

    end if 
    do inu  = 2, ace_numax
       if (l_ace_order(inu)) then 
         dim_baseB = size(base_cnlm(inu)%dicoB) * ace_kmax(inu)

         if (debug_time) time(5) = MY_MPI_WTIME()

         
         call  compute_ylm_for_neigh_ja(inu) 
         if (debug_time) then
            time_aa = MY_MPI_WTIME()
            time_ylm = time_ylm + time_aa - time(5)
         end if
         !call  compute_rad_for_neigh_ja(inu)
          if (debug_time) then
              time_bb = MY_MPI_WTIME()
              time_rad = time_rad + time_bb - time_aa
          end if
         call  compute_uniqueA_neigh_ja(inu)
          if (debug_time) then
              time_aa = MY_MPI_WTIME()
              time_unique = time_unique + time_aa - time_bb
          end if
         
         if (debug_time) then
            time(6) = MY_MPI_WTIME()
            tot_time(3) = tot_time(3) + time(6) - time(5)
         end if       
         !here should be a choice depending on body order 2, 3 - standard 
         !larger body order - adjoint, probably.    
         !call  compute_B_basis_atom(inu, baseB, d_baseB, desc_forces_local)
         call  compute_B_basis_atom_adjoint(inu, baseB, d_baseB, desc_forces_local)
         
          if (debug_time) then
            time(7) = MY_MPI_WTIME()
            tot_time(4) = tot_time(4) + time(7) - time(6)
          end if 
         !ACE_CHEM
         config_desc(iconf)%energy(icount+1+ja_pos_in_desc : icount+dim_baseB+ja_pos_in_desc, ja) = baseB(1:dim_baseB)
         !config_desc(iconf)%energy(icount+1:icount+dim_baseB, ja) = baseB(1:dim_baseB)
         if (desc_forces_local) then 
           do ix = 1,3
            !ACE_CHEM 
            config_desc(iconf)%force(icount+1+ja_pos_in_desc : icount+dim_baseB+ja_pos_in_desc, ja, 1:max_neigh_local, ix) = d_baseB(1:dim_baseB, ix, 1:max_neigh_local)
            !config_desc(iconf)%force(icount+1:icount+dim_baseB, ja, 1:max_neigh_local, ix) = d_baseB(1:dim_baseB, ix, 1:max_neigh_local)
           end do        
         end if 
         icount = icount + dim_baseB

         do iz = 2, zetaace_order
            do ii = 1, dim_baseB
              icount = icount  + 1
              !$! temp_dja = ( delta_zetaace(iz,inu) * baseB(ii)  )**(iz-1)
              !$! config_desc(iconf)%energy(icount, ja) =  temp_dja *  delta_zetaace(iz,inu) * baseB(ii)
              temp_dja = delta_zetaace(iz,inu) * ( baseB(ii)  )**(iz-1)
              !ACE_CHEM
              config_desc(iconf)%energy(icount+ja_pos_in_desc, ja) =  temp_dja *  baseB(ii)
              !config_desc(iconf)%energy(icount, ja) =  temp_dja *  baseB(ii)
            
              if (desc_forces_local) then 
                do ix = 1, 3
                  !$! config_desc(iconf)%force(icount, ja, 1:max_neigh_local, ix) =  dble(iz) * delta_zetaace(iz,inu) * &
                  !ACE_CHEM 
                  config_desc(iconf)%force(icount+ja_pos_in_desc, ja, 1:max_neigh_local, ix) =  dble(iz) * d_baseB(ii, ix, 1:max_neigh_local) * temp_dja 
                  !config_desc(iconf)%force(icount, ja, 1:max_neigh_local, ix) =  dble(iz) * d_baseB(ii, ix, 1:max_neigh_local) * temp_dja 
                end do
              end if  
            end do          
         end do


         deallocate(baseB)
         deallocate(d_baseB)
          if (debug_time) then
            time(8) = MY_MPI_WTIME()
            tot_time(5) = tot_time(5) + time(8) - time(7)
          end if 

       end if   
    end do

    if (desc_forces_local) then
      do ia_n = 1, max_neigh_local
        config_desc(iconf)%force(:, ja, 0, :) = config_desc(iconf)%force(:, ja, 0, :) - config_desc(iconf)%force(:, ja, ia_n, :)
      end do
    end if

  end do                  ! ja

  if (debug_time) then
    time(2) = MY_MPI_WTIME()
    tot_time(1) = tot_time(1) + time(2) - time(1)
  end if


  if (debug_time) then
    call repport_time(2, 0.d0, tot_time(1), "ML: full")
    call repport_time(2, 0.d0, tot_time(2), "ML: main loop: geom nn  ")
    call repport_time(2, 0.d0, tot_time(3), "ML: main loop: pre Ylm ")
    call repport_time(4, 0.d0, time_ylm, "ML: main loop: Ylm ")
    call repport_time(4, 0.d0, time_rad, "ML: main loop: radial ")
    call repport_time(4, 0.d0, time_unique, "ML: main loop: unique ") 
    call repport_time(2, 0.d0, tot_time(5), "ML: main loop: pack B ")
    call repport_time(2, 0.d0, tot_time(4), "ML: main loop: B basis ")
    call repport_time(4, 0.d0, time_amat  , "ML: inner Amat ")
    call repport_time(4, 0.d0, time_amat_00  , "ML: inner Amat 00 ")
    call repport_time(4, 0.d0, time_amat_01  , "ML: inner Amat 01 ")
    call repport_time(4, 0.d0, time_amat_02  , "ML: inner Amat 02 ")
    call repport_time(4, 0.d0, time_amat_03  , "ML: inner Amat 03 ")
    call repport_time(4, 0.d0, time_cgord  , "ML: inner cgord 00 ")

    !call repport_time(4, 0.d0, tot_time(5), "ML: nn loop: nn search")
    !call repport_time(4, 0.d0, tot_time(6), "ML: nn loop: cmm/dcmm")
    !call repport_time(6, 0.d0, tot_time(8), "ML: spher 3D : radial")
    !call repport_time(6, 0.d0, tot_time(9), "ML: spher 3D :  spher")
    !call repport_time(6, 0.d0, tot_time(10), "ML: spher 3D :    cmm")
  end if

  !stop 'end of compute_ace -----stop debug time----------'
  deallocate (xpnp)
 
 _MLD_END_ 

end subroutine compute_ace

end module module_compute_ace 


subroutine init_ace
  use iso_fortran_env,  dp => real64
  use ml_in_ndm_module, only: imm_neigh
  use module_so3, only:  radial_ace_bessel
  use module_ace_desc, only: ace_dim,  ace_lambda, ace_nmax, ace_lmax, ace_kmax, & 
                             ace_rcut_in, ace_rcut_out, ace_rcut_width_in, ace_rcut_width_out,  base_params, &
                             l_ace_order, ace_numax, delta_ace, NU_LIMIT_MAX, acenmax, acelmax, acekmax, & 
                             time_amat, time_amat_01, time_amat_02, time_amat_03, time_cgord, &
                              time_ylm, time_rad, time_unique,  &
                             time_amat_00, ace_radial_poly, l_ace_set_rcut, ace_rcut_in_list, ace_rcut_out_list, &
                             ace_rcut_width_in_list, ace_rcut_width_out_list, nkmax_order1, nkmin_order1, & 
                             zetaace_order, delta_zetaace, dim_delta_zetaace, ace_gencg, ace_chem, pos_ace_chem, & 
                             ACE_CHEM_INCOMPLETE, ACE_CHEM_STANDARD, ACE_CHEM_TS, ace_radial_chem, ace_npoints_spline
  use module_neigh_local, only: r_cut, r_cut_width,  r_cut_in, r_cut_width_in
  use module_base_cnlm, only: type_base_cnlm, base_cnlm
  use module_chemical_species, only: fix_no_of_elements
  use module_ace_radial, only: radialace, radial_spline, ACE_CHEM_RADIAL_BLOCK_HSVD, ACE_CHEM_RADIAL_HSVD, &
                               ACE_CHEM_RADIAL_RALF, ACE_CHEM_RADIAL_RANDPROJ, ACE_CHEM_RADIAL_CHEMMAP_HSVD

  use mld_logger
  use mld_mpi, only: mld_mpi_abort, mld_rank 
  implicit none
  integer :: ii, ntest, nitems2, nunit, iz, zeta_tmp 
  logical, dimension(:), allocatable  :: ltmp 

  _NAMECURRENT_("init_ace")
  _MLD_BEGIN_

  time_amat = 0.0_dp
  time_amat_00 = 0.0_dp
  time_amat_01 = 0.0_dp
  time_amat_02 = 0.0_dp
  time_amat_03 = 0.0_dp
  time_cgord = 0.0_dp
  time_ylm = 0.0_dp
  time_rad = 0.0_dp
  time_unique = 0.0_dp

  !this is delta(mui,muj) function. Probably better name will be delta_mu_ace  
  if (allocated(delta_ace)) deallocate(delta_ace) ; allocate(delta_ace(fix_no_of_elements, fix_no_of_elements))
  if ((ace_radial_chem == ACE_CHEM_RADIAL_HSVD) .or. (ace_radial_chem == ACE_CHEM_RADIAL_BLOCK_HSVD) &
      .or. (ace_radial_chem == ACE_CHEM_RADIAL_RANDPROJ) &
      .or. (ace_radial_chem == ACE_CHEM_RADIAL_CHEMMAP_HSVD)) then
    delta_ace(:,:) = 1.0_dp 
    ! acekmax = 3
  else if (ace_radial_chem == ACE_CHEM_RADIAL_RALF) then
    delta_ace(:,:) = 0.0_dp 
    do ii = 1, fix_no_of_elements
      delta_ace(ii, ii) = 1.0_dp  
    end do 
    acekmax = 1
    ace_kmax(:) = 1
  else
    call log_warning("ML: in ACE the type_chem_radial is not yet implemented.")
    stop 'stop here in init_ace'
  end if
  
  if (ace_numax> NU_LIMIT_MAX) then
    call log_critical("ML: Fatal error in "//NAMECURRENT//" concerning the number of body terms. ")
    call log_critical("The number of body terms is larger than "//vtoa(NU_LIMIT_MAX)//" , ace_numax= "//vtoa(ace_numax))
    call log_critical("possible solutions: decrese ace_numax") 
    stop 'increase ace_numax  for that descriptor. Change the sources. Ask master.'
  end if
  
  if (ace_numax<= NU_LIMIT_MAX) then 
    if (allocated(ltmp)) deallocate(ltmp); allocate(ltmp(ace_numax))
    ltmp(1:ace_numax) = l_ace_order(1:ace_numax)
    if (allocated(l_ace_order)) deallocate(l_ace_order); allocate(l_ace_order(ace_numax))
    l_ace_order(:) = ltmp(:)
    deallocate(ltmp)
  else 
    call log_critical("ML: Fatal error in "//NAMECURRENT//" concerning the number of body terms. ")
    call log_critical("The number of body terms is larger than "//vtoa(NU_LIMIT_MAX)//" , ace_numax= "//vtoa(ace_numax))
    call log_critical("possible solutions: decrese ace_numax") 
    call mld_mpi_abort("stop in "//NAMECURRENT//" increase ace_numax for ACE. Change the sources. Ask master.")
  end if   



  if (l_ace_set_rcut) then 
     if (allocated(ace_rcut_in)) deallocate(ace_rcut_in); allocate(ace_rcut_in(ace_numax))
     ntest = nitems2(ace_rcut_in_list)
     if (ntest .lt. ace_numax) then 
       call log_critical("init_ace: the number of items in ace_rcut_in_list is not equal to ace_numax")
       call log_critical("--- change and relaunch  ---")
       call mld_mpi_abort("stop in "//NAMECURRENT//" ace_rcut_in_list problem")
     end if
     read(ace_rcut_in_list, *) ace_rcut_in

      if (allocated(ace_rcut_out)) deallocate(ace_rcut_out); allocate(ace_rcut_out(ace_numax))
      ntest = nitems2(ace_rcut_out_list)
      if (ntest .lt. ace_numax) then 
        call log_critical("init_ace: the number of items in ace_rcut_out_list is not equal to ace_numax")
        call log_critical("--- change and relaunch  ---")
        call mld_mpi_abort("stop in "//NAMECURRENT//" ace_rcut_out_list problem")
      end if
      read(ace_rcut_out_list, *) ace_rcut_out

      if (allocated(ace_rcut_width_in)) deallocate(ace_rcut_width_in); allocate(ace_rcut_width_in(ace_numax))
      ntest = nitems2(ace_rcut_width_in_list)
      if (ntest .lt. ace_numax) then 
        call log_critical("init_ace: the number of items in ace_rcut_width_in_list is not equal to ace_numax")
        call log_critical("--- change and relaunch  ---")
        call mld_mpi_abort("stop in "//NAMECURRENT//" ace_rcut_width_in_list problem")
      end if
      read(ace_rcut_width_in_list, *) ace_rcut_width_in

      if (allocated(ace_rcut_width_out)) deallocate(ace_rcut_width_out); allocate(ace_rcut_width_out(ace_numax))
      ntest = nitems2(ace_rcut_width_out_list)
      if (ntest .lt. ace_numax) then 
        call log_critical("init_ace: the number of items in ace_rcut_width_out_list is not equal to ace_numax")
        call log_critical("--- change and relaunch  ---")
        call mld_mpi_abort("stop in "//NAMECURRENT//" ace_rcut_width_out_list problem")
      end if
      read(ace_rcut_width_out_list, *) ace_rcut_width_out
  else 
    if (allocated(ace_rcut_in)) deallocate(ace_rcut_in); allocate(ace_rcut_in(ace_numax))
    ace_rcut_in(:) = r_cut_in
    if (allocated(ace_rcut_out)) deallocate(ace_rcut_out); allocate(ace_rcut_out(ace_numax))
    ace_rcut_out(:) = r_cut
    if (allocated(ace_rcut_width_in)) deallocate(ace_rcut_width_in); allocate(ace_rcut_width_in(ace_numax))
    ace_rcut_width_in(:) = r_cut_width_in
    if (allocated(ace_rcut_width_out)) deallocate(ace_rcut_width_out); allocate(ace_rcut_width_out(ace_numax))
    ace_rcut_width_out(:) = r_cut_width
  end if 


  if (allocated(base_params)) deallocate(base_params); allocate(base_params(ace_numax))
  do ii = 1, ace_numax
    if (l_ace_order(ii)) then
      base_params(ii)%active = .true.  
      base_params(ii)%lmax = ace_lmax(ii)
      base_params(ii)%nmax = ace_nmax(ii)
      base_params(ii)%kmax = ace_kmax(ii)
      ! For HSVD we still enumerate chemical species, not HSVD channels.
      base_params(ii)%mumax = fix_no_of_elements
      base_params(ii)%r_cut_in = ace_rcut_in(ii)
      base_params(ii)%r_cut_out = ace_rcut_out(ii)
      base_params(ii)%r_cut_width_in = ace_rcut_width_in(ii)
      base_params(ii)%r_cut_width_out = ace_rcut_width_out(ii)
      base_params(ii)%lambda = ace_lambda(ii)
    else 
      base_params(ii)%active = .false. 
      base_params(ii)%lmax =  -1 
      base_params(ii)%nmax =  -1 
      base_params(ii)%kmax = -1
      base_params(ii)%mumax = fix_no_of_elements
      base_params(ii)%r_cut_in = -1.d0
      base_params(ii)%r_cut_out = -1.d0
      base_params(ii)%r_cut_width_in = -1.d0
      base_params(ii)%r_cut_width_in = -1.d0
      base_params(ii)%lambda = -1.d0
    end if   
  end do 


  !TODOace: a radial function for each body term. Now is unique ... 
  acenmax = maxval(base_params(:)%nmax, mask=base_params(:)%active)
  acelmax = maxval(base_params(:)%lmax, mask=base_params(:)%active)
  acekmax = maxval(base_params(:)%kmax, mask=base_params(:)%active)
  !call rfunc_ace%init(radial_pow_so3, r_cut_in, r_cut, r_cut_width_in,  r_cut_width, ace_lambda(1), acenmax)
  !TODOace endur 
  
  !TODhea! 
  call radialace%init(ace_radial_chem, ace_radial_poly, acekmax, acenmax, acelmax, fix_no_of_elements, &
                      r_cut_in, r_cut, r_cut_width_in,  r_cut_width, &
                      ace_lambda(1), ace_npoints_spline)
  !TODhea!                    

  !when_by_body! call radialace(ibody)%init(type_chem_radial, type_f_radial, &
  !when_by_body!                            base_params(ibody)%r_cut_in,  &
  !when_by_body!                            base_params(ibody)%r_cut_out, &
  !when_by_body!                            base_params(ibody)%r_cut_width_in,  &
  !when_by_body!                            base_params(ibody)%r_cut_width_out,  &
  !when_by_body!                            base_params(ibody)%lambda) 
  call radialace%build 

  ! Low-rank chemical compression is now done on-the-fly inside radialace%build
  ! when ace_chem_low_rank == 1 and ace_radial_chem == ACE_CHEM_RADIAL_HSVD.
  ! No separate post-hoc compression call needed.

  if (allocated(base_cnlm)) deallocate(base_cnlm)
  allocate(base_cnlm(ace_numax))

  ace_dim = 0
  open(newunit=nunit, file='zeta_line_number.dat', status='unknown')
  if  (mld_rank==0) write(nunit,'("# body        zeta       dim_in       dim_out")') 
  zeta_tmp = 0 
  if (l_ace_order(1)) then
    if (ace_radial_chem == ACE_CHEM_RADIAL_HSVD .or. ace_radial_chem == ACE_CHEM_RADIAL_RANDPROJ) then
      ! nkmax_order1 = acekmax
      nkmax_order1 = base_params(1)%kmax
      nkmin_order1 = 1
    else if ((ace_radial_chem == ACE_CHEM_RADIAL_BLOCK_HSVD)) then
      ! nkmax_order1 = acekmax
      nkmax_order1 = base_params(1)%kmax
      nkmin_order1 = 1
    else if (ace_radial_chem == ACE_CHEM_RADIAL_CHEMMAP_HSVD) then
      ! nkmax_order1 = acekmax
      nkmax_order1 = base_params(1)%kmax
      nkmin_order1 = 1
    else if (ace_radial_chem == ACE_CHEM_RADIAL_RALF) then
      nkmax_order1 = base_params(1)%nmax
      nkmin_order1 = 0
    end if
    ace_dim =  base_params(1)%mumax * (nkmax_order1 + 1 - nkmin_order1)
    !OLD  ace_dim =  base_params(1)%mumax * (base_params(1)%nmax+1)

    if (base_params(1)%lmax > 0) then
      call log_warning("ML: The first body term can be only with lmax = 0. The lmax will be set to 0.")
      base_params(1)%lmax = 0
    end if
    if (mld_rank==0) then
      write(nunit,'(i6, 1x, i6, 1x, i10, 1x, i10)') 1, 1, 1,  ace_dim
      !old zeta_tmp = zeta_tmp + base_params(1)%mumax * (base_params(1)%nmax+1)
      zeta_tmp = zeta_tmp + base_params(1)%mumax * (nkmax_order1 + 1 - nkmin_order1)
      if (zetaace_order > 1) then 
         do iz = 2, zetaace_order 
            write(nunit,'(i6, 1x, i6, 1x, i10, 1x, i10)') 1, iz,  zeta_tmp + 1,  zeta_tmp + base_params(1)%mumax * (base_params(1)%nmax+1)
            !old zeta_tmp = zeta_tmp + base_params(1)%mumax * (base_params(1)%nmax+1)
            zeta_tmp = zeta_tmp + base_params(1)%mumax * (nkmax_order1 + 1 - nkmin_order1)
         end do 
      end if             
    end if 
  end if  
  !wACE: case nu = 1 are treated separatelly. 
  do ii = 2, ace_numax
    if (l_ace_order(ii)) then
        if ((ace_radial_chem == ACE_CHEM_RADIAL_HSVD) .or. (ace_radial_chem == ACE_CHEM_RADIAL_BLOCK_HSVD) .or. &
          (ace_radial_chem == ACE_CHEM_RADIAL_RANDPROJ) .or. (ace_radial_chem == ACE_CHEM_RADIAL_CHEMMAP_HSVD)) then
        call base_cnlm(ii)%init(ii, 1, base_params(ii)%lmax, base_params(ii)%mumax)
      else if (ace_radial_chem == ACE_CHEM_RADIAL_RALF) then
        call base_cnlm(ii)%init(ii, base_params(ii)%nmax, base_params(ii)%lmax, base_params(ii)%mumax)
      end if
 
      !---- this is initegrated now in init ---!    
      !!build the nblb without restriction RI basis 
      !call base_cnlm(ii)%build_nblb
      !!build the nblb with permutation restriction RPI basis
      !call base_cnlm(ii)%build_nblb_same_sigma
      !----------------------------------------!

      ! start to compute the generalizedCG coefficients
      ! a. init the vectors of ordinaryCG coefficients base_cnlm(ii)%precg(j1,m1,j2,m2,j12,m12)
      call log_info("ML: Precomputing the CG coefficients with llmax "//vtoa(base_cnlm(ii)%llmax)//" for nu = "//vtoa(ii))
      call base_cnlm(ii)%build_precg(base_cnlm(ii)%llmax)
      ! b. init the vectors of RI  generalizedCG coefficients base_cnlm(ii)%dico_lb0(il)%cg(#M0)
      call log_info("ML: Precomputing the RI CG coefficients with llmax "//vtoa(base_cnlm(ii)%llmax)//" for nu = "//vtoa(ii)) 
      call base_cnlm(ii)%compute_cg(ace_gencg, ii)
      !d. build the the full index of B basis: 
      call base_cnlm(ii)%build_Bfull
      !debug call base_cnlm(ii)%dump_dicoB("dump_dicoB_fortran_"//trim(vtoa(ii))//".txt")
      call log_info("ML: ACE body "//vtoa(ii)//" has the dim_B_basis "//vtoa(size(base_cnlm(ii)%dicoB)))
      ace_dim = ace_dim + size(base_cnlm(ii)%dicoB) * base_params(ii)%kmax
      call base_cnlm(ii)%check_uniqueA_tuples
      !debug call base_cnlm(ii)%dump_uniqueA("dump_uniqueA_fortran_"//trim(vtoa(ii))//".txt")
      !debug stop 'stop here dump uniqueA'
      if (mld_rank==0) then
        write(nunit,'(i6, 1x, i6, 1x, i10, 1x, i10)') ii, 1,  zeta_tmp + 1, zeta_tmp + size(base_cnlm(ii)%dicoB) 
        zeta_tmp = zeta_tmp + size(base_cnlm(ii)%dicoB)
        if (zetaace_order > 1) then 
          do iz = 2, zetaace_order 
            write(nunit,'(i6, 1x, i6, 1x, i10, 1x, i10)') ii, iz, zeta_tmp  + 1, zeta_tmp + size(base_cnlm(ii)%dicoB)  
            zeta_tmp = zeta_tmp + size(base_cnlm(ii)%dicoB)
          end do 
        end if             
      end if 
    end if
  end do
  if (mld_rank==0) close(nunit)
    
  call print_ace_descriptor

  call log_info("ML: In this ACE descriptor the total number of B basis is "//vtoa(ace_dim))    

  if (zetaace_order > 1 ) then
    call log_info("ML: The ACE descriptor has the zeta-order ZetaACE term "//vtoa(zetaace_order)) 
    ace_dim = ace_dim * zetaace_order
  end if 

  if (zetaace_order > 1 ) then
    if (allocated(delta_zetaace)) deallocate(delta_zetaace) ; allocate(delta_zetaace(zetaace_order, ace_numax))
    ! delta_zetaace is indexed by the ZETA order (ii = 1..zetaace_order), which is INDEPENDENT
    ! of the body-order activation l_ace_order. The previous guard "if (l_ace_order(ii))" reused
    ! the body-order flag for the zeta index, so whenever zetaace_order exceeded the number of
    ! active body orders the corresponding delta_zetaace rows were left UNINITIALISED: garbage
    ! values were used in the fit while 1.0 was written to the exported LAMMPS XML, producing a
    ! fit/inference mismatch on those zeta blocks.
    !
    ! Parse EVERY zeta row (ii = 1..zetaace_order) from its own dim_delta_zetaace(ii) input.
    ! The 1.0 initialisation is only a safety net against garbage; we do NOT silently substitute
    ! 1.0 for a requested-but-missing row (that could overwrite intended/optimal non-1.0 values),
    ! we abort instead so the user supplies them explicitly.
    delta_zetaace = 1.d0
    do ii = 1, zetaace_order
      ntest = nitems2(dim_delta_zetaace(ii))
      if (ntest .lt. ace_numax) then
        call log_critical("init_ace: dim_delta_zetaace(zeta="//vtoa(ii)//") has fewer items than ace_numax = "//vtoa(ace_numax))
        call log_critical("--- provide ace_numax values for every zeta order up to zetaace_order, then relaunch ---")
        call mld_mpi_abort("stop in "//NAMECURRENT//" dim_delta_zetaace problem")
      end if
      read(dim_delta_zetaace(ii), *) delta_zetaace(ii, :)
    end do
  end if
  !r_cut_in, r_cut_out, delta_cut, lambda, ace_params
  if (r_cut_width < 0.d0 ) then
     call log_critical("r_cut_width is negative ... the mld will stop here ...")
     call mld_mpi_abort("stop in "//NAMECURRENT) 
  end if 
  if (dabs(r_cut_width_in - r_cut_width) < 1.d-16)  then
    call log_warning("r_cut_width in and out have different values ...") 
    call log_warning("Both should be equal in ACE formalism. The \in\ value will be taken the same as \out\ value. ")
    r_cut_width_in = r_cut_width
  end if



  select case (ace_chem)
     case(ACE_CHEM_INCOMPLETE)
       pos_ace_chem = 0 
     case (ACE_CHEM_STANDARD)
       pos_ace_chem = ace_dim 
     case (ACE_CHEM_TS)  
       call log_critical("ML: ACE_CHEM_TS ace_chem=2 is not implemented yet. ")
       call log_critical("ML: Please, use ACE_CHEM_INCOMPLETE (ace_chem=0) or ACE_CHEM_COMPLETE (ace_chem=1) ")
       call mld_mpi_abort("stop in "//NAMECURRENT//" ACE_CHEM_TS is not implemented yet. ")
     case default
       call log_critical("ML: ACE_CHEM is not recognized. ace_chem can be 0,1 or 2. ")
       call mld_mpi_abort("stop in "//NAMECURRENT//" ACE_CHEM is not recognized. ")
  end select

  _MLD_END_ 
end subroutine init_ace
