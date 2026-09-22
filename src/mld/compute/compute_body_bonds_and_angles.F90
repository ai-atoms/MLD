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

module module_compute_body_bonds_and_angles

contains 
subroutine compute_body_bonds_and_angles(i_start_at, i_final_at, d_n_neigh, d_kind_neigh, iconf)

  USE module_kind_variables, ONLY: double
#ifdef MLD_NDM
  use gen_com_m, ONLY: A2cm, lperiod
  use gen_com_m_ml, ONLY: imm, bg, at
  use tab_imm_m_ml, ONLY: xp, iwmax2 ! xp utilise ??
#else
  use ondm_gen_com_m, ONLY: imm, lperiod
  use ondm_tab_imm_m, ONLY: iwmax2, xp
#endif
  use ml_in_ndm_module, ONLY: imm_neigh, desc_forces, rangml
  use module_body_desc, only: l_body_order, desc_forces_bond, &
                              dim_desc_body, time_body, time_invariants, time_polynomial, &
                              tmp_force, tmp_real
  use module_neigh_local, only: r_cut, ja_atom, iw2, max_neigh_local, i_type, i_type_db, i_central, r_central,  &
                                ur_central,  tmp_dxp, tmp_xp, reallocate_neigh_ja, build_local_neighbours_ja
#ifdef MLD_NDM                               
  use notperiod_mod
#endif


  use derived_types, only: config_real, config_desc
  use time_check_general, only: time, debug_time, MY_MPI_WTIME
#ifndef MLD_NDM
  use ondm_transform_coord, only: ondm_notperiod
#endif

  implicit none

  integer, intent(in)  :: i_start_at, i_final_at
  integer, dimension(imm), intent(out)   :: d_n_neigh
  integer, dimension(imm, imm_neigh), intent(out)    :: d_kind_neigh
  ! real(double), dimension(:), allocatable ::   tmp_real
  ! double precision,dimension(dim_xdesc, imm),intent(out)   :: config_desc(iconf)%energy
  ! double precision,dimension(dim_xdesc, imm,0:imm_neigh, 3),intent(out) :: config_desc(iconf)%force
  integer, optional    :: iconf

  logical  :: small
  real(double), dimension(:, :), allocatable   :: xpnp
  real(double) :: t_nn,  t_allo 
  integer  :: ia, ja, ix
  ! double precision :: factor_ia, factor_ja
  integer  :: max_neigh, icnt, ii
  logical  :: desc_forces_local


  max_neigh = imm_neigh
  time_polynomial(:) = 0.d0
  time_invariants(:) = 0.d0

  if ((i_start_at == 0) .and. (i_final_at == 0)) then
    d_n_neigh(:) = 0
    d_kind_neigh(:, :) = 0
    config_desc(iconf)%energy(:, :) = 0.d0
    ! config_desc(iconf)%force(:,:,:,:)=0.d0
    return
  end if
  desc_forces_local = desc_forces .and. (config_real(iconf)%has_force .or. config_real(iconf)%has_stress)
  desc_forces_bond = desc_forces_local
  small = .false.
  if (present(iconf)) then
    small = config_real(iconf)%small
  end if
!xpnp inutile  
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
  ! call cryst_to_cart (imm, xpnp, bg, -1)
  t_allo = 0.d0
  t_nn = 0.d0
  d_n_neigh(:) = 0
  d_kind_neigh(:, :) = 0
  config_desc(iconf)%energy(:, :) = 0.d0
  if (desc_forces_local) config_desc(iconf)%force(:, :, :, :) = 0.d0
  if (i_start_at == 1) iw2 = 0
  if (i_start_at > 1) iw2 = iwmax2(i_start_at - 1)

  do ja = i_start_at, i_final_at

    if (debug_time) time(1) = MY_MPI_WTIME()
    ! provide r_central, i_central, i_type, tmp_dxp, tmp_xp, d_n_neigh, max_neigh_local
    ja_atom = ja
    !call local_neighbours(          iconf, ja, imm, xpnp, r_cut, d_n_neigh, d_kind_neigh)
    call build_local_neighbours_ja( iconf, ja, imm, xpnp, r_cut, d_n_neigh, d_kind_neigh, &
                                    r_central, i_type, i_type_db, i_central, tmp_dxp, tmp_xp, &
                                    max_neigh_local, iw2)

    if (max_neigh_local == 0) cycle
    if (max_neigh_local > max_neigh) then
      if (rangml == 0) then
        write (6, *) 'ML: Fatal error in compute_cluster_bonds_and_angles the number of atoms. '
        write (6, *) 'The number of max_neig_local is bigger than max_neigh r_cut, max_neigh_local, max_neigh', r_cut, max_neigh, max_neigh_local
        write (6, *) 'possible solutions: decrease r_cut or increase max_neigh'
      end if
      stop 'decrease rcut for that descriptor'
    end if

    if (debug_time) then
      time(2) = MY_MPI_WTIME()
      t_nn = t_nn + time(2) - time(1)
    end if

    !call reallocate_neigh_ja()
    call reallocate_neigh_ja (max_neigh_local, ja_atom, i_central, i_type, i_type_db, r_central, tmp_dxp)
    call compute_cos_and_fcut_ja(desc_forces_local)
    call compute_ur_transformed_distances(desc_forces_local)

    icnt = 0
    time(3) = MY_MPI_WTIME()
    t_allo = t_allo + time(3) - time(2)
    if (l_body_order(2)) then
      if (allocated(tmp_real)) deallocate (tmp_real); allocate (tmp_real(dim_desc_body(2)))
      if (desc_forces_local) then
        if (allocated(tmp_force)) deallocate (tmp_force); allocate (tmp_force(dim_desc_body(2), max_neigh_local, 3))
      end if
      call bond_order_2()
      do ii = 1, dim_desc_body(2)
        icnt = icnt + 1
        config_desc(iconf)%energy(icnt, ja) = tmp_real(ii)
        if (desc_forces_local) then
          do ix = 1, 3
            config_desc(iconf)%force(icnt, ja, 1:max_neigh_local, ix) = tmp_force(ii, 1:max_neigh_local, ix)
            ! write (*,*) SUM(config_desc(iconf)%force(icnt,ja,:,ix))
          end do
        end if
      end do
    end if

    if (debug_time) then
      time(4) = MY_MPI_WTIME()
      time_body(2) = time_body(2) + time(4) - time(3)
    end if

    if (l_body_order(3)) then
      if (allocated(tmp_real)) deallocate (tmp_real); allocate (tmp_real(dim_desc_body(3)))
      if (desc_forces_local) then
        if (allocated(tmp_force)) deallocate (tmp_force); allocate (tmp_force(dim_desc_body(3), max_neigh_local, 3))
      end if
      call bond_order_3()
      do ii = 1, dim_desc_body(3)
        icnt = icnt + 1
        config_desc(iconf)%energy(icnt, ja) = tmp_real(ii)
        if (desc_forces_local) then
          do ix = 1, 3
            config_desc(iconf)%force(icnt, ja, 1:max_neigh_local, ix) = tmp_force(ii, 1:max_neigh_local, ix)
            ! write (*,*) SUM(config_desc(iconf)%force(icnt,ja,:,ix))
          end do
        end if
      end do
    end if

    if (debug_time) then
      time(5) = MY_MPI_WTIME()
      time_body(3) = time_body(3) + time(5) - time(4)
    end if

    if (l_body_order(4)) then
      if (allocated(tmp_real)) deallocate (tmp_real); allocate (tmp_real(dim_desc_body(4)))
      if (desc_forces_local) then
        if (allocated(tmp_force)) deallocate (tmp_force); allocate (tmp_force(dim_desc_body(4), max_neigh_local, 3))
      end if
      call bond_order_4()
      do ii = 1, dim_desc_body(4)
        icnt = icnt + 1
        config_desc(iconf)%energy(icnt, ja) = tmp_real(ii)
        if (desc_forces_local) then
          do ix = 1, 3
            config_desc(iconf)%force(icnt, ja, 1:max_neigh_local, ix) = tmp_force(ii, 1:max_neigh_local, ix)
            ! write (*,*) SUM(config_desc(iconf)%force(icnt,ja,:,ix))
          end do
        end if
      end do
    end if

    if (debug_time) then
      time(6) = MY_MPI_WTIME()
      time_body(4) = time_body(4) + time(6) - time(5)
    end if

    if (desc_forces_local) then
      do ia = 1, max_neigh_local
        config_desc(iconf)%force(:, ja, 0, :) = config_desc(iconf)%force(:, ja, 0, :) - config_desc(iconf)%force(:, ja, ia, :)
      end do
    end if

  end do                  ! ja main loop

  if (debug_time) then
    call repport_time(2, 0.d0, t_nn, "ML: time for search of nn")
    call repport_time(2, 0.d0, t_allo, "ML: time for allocation post nn ")
    call repport_time(2, 0.d0, time_body(2), "ML: time for order 2")
    call repport_time(2, 0.d0, time_body(3), "ML: time for order 3")
    call repport_time(2, 0.d0, time_body(4), "ML: time for order 4")
    call repport_time(4, 0.d0, time_invariants(4), "ML: time 4 invariants ")
    call repport_time(4, 0.d0, time_polynomial(4), "ML: time 4 polynoms  ")
  end if
  deallocate (xpnp, ur_central, r_central, i_central, i_type, tmp_real, tmp_dxp)

end subroutine compute_body_bonds_and_angles

end module module_compute_body_bonds_and_angles

subroutine compute_ur_transformed_distances(desc_forces_local)
  use module_neigh_local, only: max_neigh_local, ur_central, d_ur_central, r_central, &
                                tmpcos_dxp
  use module_body_desc, only: bond_dist_pure, bond_dist_exp, bond_dist_inverse, bond_dist_transform, &
                              bond_beta, bond_dist_ann
  implicit none

  logical  :: desc_forces_local
  integer  :: ii

  ! compute ur_central ... the transformation of r_central ...
  if (allocated(ur_central)) deallocate (ur_central); allocate (ur_central(max_neigh_local))
  ! transform the radial function ...
  if (desc_forces_local) then
    if (allocated(d_ur_central)) deallocate (d_ur_central); allocate (d_ur_central(3, max_neigh_local))
  end if
  select case (bond_dist_transform)
  case (bond_dist_pure)
    ur_central(1:max_neigh_local) = r_central(1:max_neigh_local)
    if (desc_forces_local) then
      do ii = 1, 3
        d_ur_central(ii, 1:max_neigh_local) = tmpcos_dxp(ii, 1:max_neigh_local)
      end do
    end if
  case (bond_dist_exp)
    ur_central(1:max_neigh_local) = dexp(-bond_beta*r_central(1:max_neigh_local))
    if (desc_forces_local) then
      do ii = 1, 3
        d_ur_central(ii, 1:max_neigh_local) = tmpcos_dxp(ii, 1:max_neigh_local)*(-bond_beta)*ur_central(1:max_neigh_local)
      end do
    end if
  case (bond_dist_inverse)
    ur_central(1:max_neigh_local) = bond_dist_ann**bond_beta/r_central(1:max_neigh_local)**bond_beta
    if (desc_forces_local) then
      do ii = 1, 3
        d_ur_central(ii, 1:max_neigh_local) = (-bond_beta)*tmpcos_dxp(ii, 1:max_neigh_local)*ur_central(1:max_neigh_local)/r_central(1:max_neigh_local)
      end do
    end if
  end select

end subroutine compute_ur_transformed_distances



subroutine bond_order_2()

  use module_neigh_local, only: max_neigh_local, &
                                ur_central, d_ur_central, &
                                r_fcut, d_r_fcut
  use module_body_desc, only: basis2, dim_desc_body, &
                              tmp_force, tmp_real, desc_forces_bond
  implicit none

  integer  :: i2, id
  real(kind(1.d0)), dimension(1)   :: adp2, ads2
  real(kind(1.d0)), dimension(3, 1)      :: d_adp2, d_ads2

  real(kind(1.d0))     :: tmp_cut, tmp_pol, x1
  real(kind(1.d0)), dimension(3)   :: d_tmp_cut, d_tmp_pol
  real(kind(1.d0)), dimension(3, dim_desc_body(2))   :: d_tmp_real
  integer  :: ix


  ! compute all the body 2 distance
  ! j--x1--i2 2
  if (max_neigh_local <= 1) return
  tmp_real(:) = 0.d0
  if (desc_forces_bond) then
    d_tmp_real(:, :) = 0.d0
    tmp_force(:, :, :) = 0.d0
  end if

  do i2 = 1, max_neigh_local
    x1 = ur_central(i2)
    adp2(1) = x1
    ads2(1) = 1.d0

    if (desc_forces_bond) then
      d_adp2(1:3, 1) = d_ur_central(1:3, i2)
      d_ads2(1:3, 1) = 0.d0
      d_tmp_cut(1:3) = d_r_fcut(1:3, i2)
    end if
    tmp_cut = r_fcut(i2)

    !an if here ... for PiP
    
    !and else for nbody kernels ... 

    do id = 1, dim_desc_body(2)
      tmp_pol = adp2(1)**basis2(id)%ka(1)
      tmp_real(id) = tmp_real(id) + tmp_cut*ads2(basis2(id)%b)*tmp_pol
      if (desc_forces_bond) then
        if (basis2(id)%ka(1) == 0) then
          d_tmp_pol(1:3) = 0
        else
          d_tmp_pol(1:3) = basis2(id)%ka(1)*tmp_pol/adp2(1)*d_adp2(1:3, 1)
        end if
        ! d_tmp_real(1:3,id) =  d_tmp_real(1:3,id) +  ads2(basis2(id)%b)*(d_tmp_cut(1:3)*tmp_pol + tmp_cut*d_tmp_pol(1:3))
        d_tmp_real(1:3, id) = ads2(basis2(id)%b)*(d_tmp_cut(1:3)*tmp_pol + tmp_cut*d_tmp_pol(1:3))
      end if
    end do

    if (desc_forces_bond) then
      do ix = 1, 3
        tmp_force(1:dim_desc_body(2), i2, ix) = d_tmp_real(ix, 1:dim_desc_body(2))/2.d0
      end do
    end if
  end do

  tmp_real(:) = tmp_real(:)/2.d0                   ! because i1 there are no symmetries ...
  ! write (*,*) d_tmp_real(1,:)

end subroutine bond_order_2



subroutine bond_order_3()
  ! compute all the body3 triangle with the following structure.
  ! Please pay attention that the distances are already r_cut -ed.
  ! j--x1--i2 2
  ! \ x3
  !   x2
  !    \
  !     i3   3

  use module_neigh_local, only: max_neigh_local, i_type, &
                                r_central, &
                                ur_central, d_ur_central, &
                                r_fcut, d_r_fcut, &
                                tmp_dxp
  
  use module_body_desc, only: basis3, dim_desc_body, &
                              tmp_force, tmp_real, desc_forces_bond
  implicit none

  integer  :: i2, i3, id

  real(kind(1.d0)), dimension(3)   :: adp3
  real(kind(1.d0)), dimension(3, 3)      :: d2_adp3, d3_adp3
  real(kind(1.d0)), dimension(1)   :: ads3
  real(kind(1.d0)), dimension(3, 1)      :: d2_ads3, d3_ads3

  real(kind(1.d0)), dimension(dim_desc_body(3), max_neigh_local, 3)    :: d2_tmp_real, d3_tmp_real
  real(kind(1.d0)), dimension(3)   :: d2_tmp_pol, d3_tmp_pol, &
                                      d2_tmp_cut, d3_tmp_cut, &
                                      d2_x1, d2_x2, d2_x3, &
                                      d3_x1, d3_x2, d3_x3
  real(kind(1.d0))     :: k1, k2, k3, t1, t2, t3, t123

  real(kind(1.d0))     :: cos_2j3, tmp_cut, tmp_pol, x1, x2, x3


  if (max_neigh_local <= 1) return
  tmp_real(:) = 0.d0
  if (desc_forces_bond) then
    d2_tmp_real(:, :, :) = 0.d0
    d3_tmp_real(:, :, :) = 0.d0
    tmp_force(:, :, :) = 0.d0
  end if
  do i2 = 1, max_neigh_local - 1
    x1 = ur_central(i2)
    if (desc_forces_bond) then
      d2_x1(1:3) = d_ur_central(1:3, i2)
      d3_x1(1:3) = 0.d0
    end if
    do i3 = i2 + 1, max_neigh_local
      x2 = ur_central(i3)
      cos_2j3 = dot_product(tmp_dxp(1:3, i2), tmp_dxp(1:3, i3))/(r_central(i2)*r_central(i3))
      x3 = cos_2j3
      if (desc_forces_bond) then
        d2_x2(1:3) = 0.d0
        d3_x2(1:3) = d_ur_central(1:3, i3)
        d2_x3(1:3) = tmp_dxp(1:3, i3)/(r_central(i2)*r_central(i3)) - tmp_dxp(1:3, i2)*cos_2j3/r_central(i2)**2
        d3_x3(1:3) = tmp_dxp(1:3, i2)/(r_central(i2)*r_central(i3)) - tmp_dxp(1:3, i3)*cos_2j3/r_central(i3)**2
      end if
      ! AA case - the same species ...
      if (i_type(i2) == i_type(i3)) then

        adp3(1) = x1 + x2
        adp3(2) = x1*x2
        adp3(3) = x3
        ads3(1) = 1.d0
        if (desc_forces_bond) then
          d2_adp3(1:3, 1) = d2_x1(1:3)
          d2_adp3(1:3, 2) = d2_x1(1:3)*x2
          d2_adp3(1:3, 3) = d2_x3(1:3)
          d3_adp3(1:3, 1) = d3_x2(1:3)
          d3_adp3(1:3, 2) = d3_x2(1:3)*x1
          d3_adp3(1:3, 3) = d3_x3(1:3)
          d2_ads3(1:3, 1) = 0.d0
          d3_ads3(1:3, 1) = 0.d0
        end if

      else

        ! AB case
        ! if (i_type(i2)/=i_type(i3)) then
        adp3(1) = x1
        adp3(2) = x2
        adp3(3) = x3
        ads3(1) = 1.d0
        if (desc_forces_bond) then
          d2_adp3(1:3, 1) = d2_x1(1:3)
          d2_adp3(1:3, 2) = d2_x2(1:3)
          d2_adp3(1:3, 3) = d2_x3(1:3)
          d2_ads3(1:3, 1) = 0.d0
          d3_adp3(1:3, 1) = d3_x1(1:3)
          d3_adp3(1:3, 2) = d3_x2(1:3)
          d3_adp3(1:3, 3) = d3_x3(1:3)
          d3_ads3(1:3, 1) = 0.d0
        end if

      end if

      tmp_cut = r_fcut(i2)*r_fcut(i3)
      if (desc_forces_bond) then
        d2_tmp_cut(1:3) = d_r_fcut(1:3, i2)*r_fcut(i3)
        d3_tmp_cut(1:3) = r_fcut(i2)*d_r_fcut(1:3, i3)
      end if

      do id = 1, dim_desc_body(3)
        k1 = dble(basis3(id)%ka(1))
        k2 = dble(basis3(id)%ka(2))
        k3 = dble(basis3(id)%ka(3))
        t1 = adp3(1)**(k1 - 1)
        t2 = adp3(2)**(k2 - 1)
        if (adp3(3) == 0.d0) then
          t3 = 0.d0
        else
          t3 = adp3(3)**(k3 - 1)
        end if
        t123 = t1*t2*t3
        ! tmp_pol= PRODUCT(adp3(:)**basis3(id)%ka(:))
        tmp_pol = t123*adp3(1)*adp3(2)*adp3(3)
        tmp_real(id) = tmp_real(id) + tmp_cut*ads3(basis3(id)%b)*tmp_pol
        if (desc_forces_bond) then
          d2_tmp_pol(1:3) = basis3(id)%ka(1)*t123*adp3(2)*adp3(3)*d2_adp3(1:3, 1) + &
                            ! this one is zero: basis3(id)%ka(3)*t123*adp3(1)*adp3(3)*d2_adp3(1:3,2) + &
                            basis3(id)%ka(3)*t123*adp3(1)*adp3(2)*d2_adp3(1:3, 3)
          d2_tmp_real(id, i2, 1:3) = d2_tmp_real(id, i2, 1:3) + d2_tmp_cut(1:3)*ads3(basis3(id)%b)*tmp_pol + &
                                     tmp_cut*ads3(basis3(id)%b)*d2_tmp_pol(1:3)



          d3_tmp_pol(1:3) = basis3(id)%ka(2)*t123*adp3(1)*adp3(3)*d3_adp3(1:3, 2) + &
                            ! this one is zero: basis3(id)%ka(1)*t123*adp3(2)*adp3(3)*d3_adp3(1:3,1) + &
                            basis3(id)%ka(3)*t123*adp3(1)*adp3(2)*d3_adp3(1:3, 3)
          d3_tmp_real(id, i3, 1:3) = d3_tmp_real(id, i3, 1:3) + d3_tmp_cut(1:3)*ads3(basis3(id)%b)*tmp_pol + &
                                     tmp_cut*ads3(basis3(id)%b)*d3_tmp_pol(1:3)

        end if
      end do                  ! dim_desc
    end do                  ! i3
  end do                  ! i2

  if (desc_forces_bond) then
    tmp_force(1:dim_desc_body(3), 1:max_neigh_local, 1:3) = (d2_tmp_real(1:dim_desc_body(3), 1:max_neigh_local, 1:3) + &
                                                             d3_tmp_real(1:dim_desc_body(3), 1:max_neigh_local, 1:3))/3.d0
  end if
  tmp_real(:) = tmp_real(:)/3.d0                   ! because only i2 and i3 are symmetric indexes and j is summed on all procs ...
end subroutine bond_order_3



subroutine bond_order_4()
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

  !use, intrinsic :: ISO_Fortran_env, only: REAL64
  use module_neigh_local, only: max_neigh_local, &
                                r_central, &
                                ur_central, d_ur_central, &
                                r_fcut, d_r_fcut, &
                                tmp_dxp
  use module_body_desc, only: dim_desc_body, &
                              tmp_force, tmp_real, &
                              body_D_max, time_invariants, time_polynomial, bb4, kk4, desc_forces_bond
  use time_check_general, only: debug_time, MY_MPI_WTIME

  implicit none

  integer  :: i2, i3, i4, id, ip, im, ix
  real(kind(1.d0)), dimension(6)   :: adp4
  real(kind(1.d0)), dimension(3, 6)      :: d2_adp4, d3_adp4, d4_adp4
  real(kind(1.d0)), dimension(12)  :: ads4
  real(kind(1.d0)), dimension(3, 12)     :: d2_ads4, d3_ads4, d4_ads4
  real(kind(1.d0)), dimension(dim_desc_body(4), max_neigh_local, 3)    :: d2_tmp_real, d3_tmp_real, d4_tmp_real
  real(kind(1.d0)), dimension(dim_desc_body(4), 3)   :: tmp_i2, tmp_i3, tmp_i4
  real(kind(1.d0))     :: time1, time2, time3
  real(kind(1.d0))     :: pmatrix(6, 0:body_D_max(4)), dd_pmatrix(6, 0:body_D_max(4)), d2_pmatrix(3, 6, 0:body_D_max(4)), d3_pmatrix(3, 6, 0:body_D_max(4)), d4_pmatrix(3, 6, 0:body_D_max(4))
  integer, dimension(:), pointer   :: iidk
  real(kind(1.d0))     :: cos_2j3, cos_3j4, cos_2j4, tmp_cut, tmp_pol, x1, x2, x3, x4, x5, x6, &
                          x1_2, x1_3, x2_2, x2_3, x3_2, x3_3, &
                          p1, p2, p3, p4, p5, p6, &
                          pm1, pm2, pm3, pm4, pm5, pm6, pm123, pm456, at, tt, ta, &
                          d1, d2, d3, d4, d5, d6

  real(kind(1.d0)), dimension(3)   :: d2_tmp_pol, d3_tmp_pol, d4_tmp_pol, &
                                      d2_tmp_cut, d3_tmp_cut, d4_tmp_cut, &
                                      d2_x1, d2_x2, d2_x3, d2_x4, d2_x5, d2_x6, &
                                      d3_x1, d3_x2, d3_x3, d3_x4, d3_x5, d3_x6, &
                                      d4_x1, d4_x2, d4_x3, d4_x4, d4_x5, d4_x6
  ! real(kind(1.d0)) :: ddot
  ! real(REAL64), external    :: ddot


  if (max_neigh_local <= 1) return
  tmp_real(:) = 0.d0
  if (desc_forces_bond) then
    d2_tmp_real(:, :, :) = 0.d0
    d3_tmp_real(:, :, :) = 0.d0
    d4_tmp_real(:, :, :) = 0.d0
    tmp_force(:, :, :) = 0.d0

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

  do i2 = 1, max_neigh_local - 2
    x1 = ur_central(i2)
    x1_2 = x1**2
    x1_3 = x1_2*x1
    if (desc_forces_bond) then
      d2_x1(1:3) = d_ur_central(1:3, i2)
    end if

    do i3 = i2 + 1, max_neigh_local - 1
      x2 = ur_central(i3)
      x2_2 = x2**2
      x2_3 = x2_2*x2
      cos_2j3 = dot_product(tmp_dxp(:, i2), tmp_dxp(:, i3))/(r_central(i2)*r_central(i3))
      x4 = cos_2j3
      if (desc_forces_bond) then
        d3_x2(1:3) = d_ur_central(1:3, i3)
        d2_x4(1:3) = tmp_dxp(1:3, i3)/(r_central(i2)*r_central(i3)) - tmp_dxp(1:3, i2)*cos_2j3/r_central(i2)**2
        d3_x4(1:3) = tmp_dxp(1:3, i2)/(r_central(i2)*r_central(i3)) - tmp_dxp(1:3, i3)*cos_2j3/r_central(i3)**2
      end if

      do i4 = i3 + 1, max_neigh_local
        x3 = ur_central(i4)
        x3_2 = x3**2
        x3_3 = x3_2*x3
        cos_3j4 = dot_product(tmp_dxp(:, i3), tmp_dxp(:, i4))/(r_central(i3)*r_central(i4))
        x6 = cos_3j4

        cos_2j4 = dot_product(tmp_dxp(:, i2), tmp_dxp(:, i4))/(r_central(i2)*r_central(i4))
        x5 = cos_2j4
        if (desc_forces_bond) then
          d4_x3(1:3) = d_ur_central(1:3, i4)
          d3_x6(1:3) = tmp_dxp(1:3, i4)/(r_central(i3)*r_central(i4)) - tmp_dxp(1:3, i3)*cos_3j4/r_central(i3)**2
          d4_x6(1:3) = tmp_dxp(1:3, i3)/(r_central(i3)*r_central(i4)) - tmp_dxp(1:3, i4)*cos_3j4/r_central(i4)**2
          d2_x5(1:3) = tmp_dxp(1:3, i4)/(r_central(i2)*r_central(i4)) - tmp_dxp(1:3, i2)*cos_2j4/r_central(i2)**2
          d4_x5(1:3) = tmp_dxp(1:3, i2)/(r_central(i2)*r_central(i4)) - tmp_dxp(1:3, i4)*cos_2j4/r_central(i4)**2
        end if
        if (debug_time) time1 = MY_MPI_WTIME()
        ! AAA case
        ! if ( (i_type(i2)==i_type(i3)).and.(i_type(i2)==i_type(i4)) ) then
        adp4(1) = x1 + x2 + x3
        adp4(2) = x4 + x5 + x6
        !adp4(3) = x1**2 + x2**2 + x3**2
        adp4(3) = x1_2 + x2_2 + x3_2
        adp4(4) = x1*x4 + x2*x5 + x3*x6
        adp4(5) = x4**3 + x5**3 + x6**3
        ! adp4(6) = x1**3 + x2**3 + x3**3 + x4**2*x5 + x4*x6**2 + x6*x5**2
        adp4(6) = x1_3 + x2_3 + x3_3 + x4**2*x5 + x4*x6**2 + x6*x5**2

        if (desc_forces_bond) then

          ! 2 is in 1 - dist ;  4,5 - angles
          d2_adp4(1:3, 1) = d2_x1(1:3)
          d2_adp4(1:3, 2) = d2_x4(1:3) + d2_x5(1:3)
          d2_adp4(1:3, 3) = 2.d0*x1*d2_x1(1:3)
          d2_adp4(1:3, 4) = d2_x1(1:3)*x4 + x1*d2_x4(:3) + x2*d2_x5(1:3)
          d2_adp4(1:3, 5) = 3.d0*x4**2*d2_x4(1:3) + 3.d0*x5**2*d2_x5(1:3)
          d2_adp4(1:3, 6) = 3.d0*x1**2*d2_x1(1:3) + 2.d0*x4*d2_x4(1:3)*x5 + x4**2*d2_x5(1:3) + d2_x4(1:3)*x6**2 + 2.d0*x5*x6*d2_x5(1:3)

          ! 3 is in 2 - dist ;  4,6 - angles
          d3_adp4(1:3, 1) = d3_x2(1:3)
          d3_adp4(1:3, 2) = d3_x4(1:3) + d3_x6(1:3)
          d3_adp4(1:3, 3) = 2.d0*x2*d3_x2(1:3)
          d3_adp4(1:3, 4) = x1*d3_x4(1:3) + d3_x2(1:3)*x5 + x3*d3_x6(1:3)
          d3_adp4(1:3, 5) = 3.d0*x4**2*d3_x4(1:3) + 3.d0*x6**2*d3_x6(1:3)
          d3_adp4(1:3, 6) = 3.d0*x2**2*d3_x2(1:3) + 2.d0*x4*d3_x4(1:3)*x5 + d3_x4(1:3)*x6**2 + x4*2.d0*x6*d3_x6(1:3) + d3_x6(1:3)*x5**2

          ! 4 is in 3 - dist ;  5,6 - angles
          d4_adp4(1:3, 1) = d4_x3(1:3)
          d4_adp4(1:3, 2) = d4_x5(1:3) + d4_x6(1:3)
          d4_adp4(1:3, 3) = 2.d0*x3*d4_x3(1:3)
          d4_adp4(1:3, 4) = x2*d4_x5(1:3) + x3*d4_x6(1:3) + d4_x3(1:3)*x6
          d4_adp4(1:3, 5) = 3.d0*x5**2*d4_x5(1:3) + 3.d0*x6**2*d4_x6(1:3)
          d4_adp4(1:3, 6) = 3.d0*x3**2*d4_x3(1:3) + x4**2*d4_x5(1:3) + 2.d0*x4*x6*d4_x6(1:3) + d4_x6(1:3)*x5**2 + 2.d0*x6*x5*d4_x5(1:3)

        end if

        ads4(1) = 1.d0
        ads4(2) = x1*x5 + x2*x4 + x3*x6
        ads4(3) = x4**2 + x5**2 + x6**2
        ads4(4) = x1_2*x3 + x1*x2_2 + x2*x3_2
        ads4(5) = x1*x2*x4 + x1*x3*x6 + x2*x3*x6
        ads4(6) = x1*x4**2 + x2*x6**2 + x3*x5**2
        ads4(7) = x1_2*x4 + x2_2*x6 + x3_2*x5
        ads4(8) = x1*x4*x5 + x2*x4*x6 + x3*x5*x6
        ads4(9) = x4**2*x6 + x4*x5**2 + x5*x6**2
        ads4(10) = ads4(2)*ads4(3)
        ads4(11) = ads4(2)**2
        ads4(12) = ads4(3)**2

        if (desc_forces_bond) then
          ! 2 is in 1 - dist ;  4,5 - angles
          d2_ads4(1:3, 1) = 0.d0
          d2_ads4(1:3, 2) = d2_x1(1:3)*x5 + x1*d2_x5(1:3) + x2*d2_x4(1:3)
          d2_ads4(1:3, 3) = 2.d0*x4*d2_x4(1:3) + 2.d0*x5*d2_x5(1:3)
          d2_ads4(1:3, 4) = 2.d0*x1*x3*d2_x1(1:3) + d2_x1(1:3)*x2_2
          d2_ads4(1:3, 5) = d2_x1(1:3)*x2*x4 + x1*x2*d2_x4(1:3) + d2_x1(1:3)*x3*x6
          d2_ads4(1:3, 6) = d2_x1(1:3)*x4**2 + x1*2.d0*x4*d2_x4(1:3) + x3_2*d2_x5(1:3)
          d2_ads4(1:3, 7) = 2.d0*x1*d2_x1(1:3)*x4 + x1**2*d2_x4(1:3) + x3_2*d2_x5(1:3)
          d2_ads4(1:3, 8) = d2_x1(1:3)*x4*x5 + x1*d2_x4(1:3)*x5 + x1*x4*d2_x5(1:3) + x2*x6*d2_x4(1:3) + x3*x6*d2_x5(1:3)
          d2_ads4(1:3, 9) = 2.d0*x4*d2_x4(1:3)*x6 + d2_x4(1:3)*x5**2 + 2.d0*x5*x4*d2_x5(1:3) + d2_x5(1:3)*x6**2
          d2_ads4(1:3, 10) = d2_ads4(1:3, 2)*ads4(3) + ads4(2)*d2_ads4(1:3, 3)
          d2_ads4(1:3, 11) = 2.d0*ads4(2)*d2_ads4(1:3, 2)
          d2_ads4(1:3, 12) = 2.d0*ads4(3)*d2_ads4(1:3, 3)

          ! 3 is in 2 - dist ;  4,6 - angles
          d3_ads4(1:3, 1) = 0.d0
          d3_ads4(1:3, 2) = d3_x2(1:3)*x4 + d3_x4(1:3)*x2 + x3*d3_x6(1:3)
          d3_ads4(1:3, 3) = 2.d0*x4*d3_x4(1:3) + 2.d0*x6*d3_x6(1:3)
          d3_ads4(1:3, 4) = 2.d0*x1*x2*d3_x1(1:3) + d3_x2(1:3)*x3_2
          d3_ads4(1:3, 5) = x1*x4*d3_x2(1:3) + x1*x2*d3_x4(1:3) + x1*x3*d3_x6(1:3) + d3_x2(1:3)*x3*x6 + x2*x3*d3_x6(1:3)
          d3_ads4(1:3, 6) = 2.d0*x1*x4*d3_x4(1:3) + 2.d0*x2*x6*d3_x6(1:3)
          d3_ads4(1:3, 7) = d3_x4(1:3)*x1_2 + 2.d0*x2*x6*d3_x2(1:3) + x2_2*d3_x6(1:3)
          d3_ads4(1:3, 8) = x1*x5*d3_x4(1:3) + d3_x2(1:3)*x4*x6 + d3_x4(1:3)*x2*x6 + d3_x6(1:3)*x4*x2 + x3*x5*d3_x6(1:3)
          d3_ads4(1:3, 9) = 2.d0*x4*x6*d3_x4(1:3) + x4**2*d3_x6(1:3) + d3_x4(1:3)*x5**2 + 2.d0*x5*x6*d3_x6(1:3)
          d3_ads4(1:3, 10) = d3_ads4(1:3, 2)*ads4(3) + ads4(2)*d3_ads4(1:3, 3)
          d3_ads4(1:3, 11) = 2.d0*ads4(2)*d3_ads4(1:3, 2)
          d3_ads4(1:3, 12) = 2.d0*ads4(3)*d3_ads4(1:3, 3)

          ! 4 is in 3 - dist ;  5,6 - angles
          d4_ads4(1:3, 1) = 0.d0
          d4_ads4(1:3, 2) = x1*d4_x5(1:3) + d4_x3(1:3)*x6 + x3*d4_x6(1:3)
          d4_ads4(1:3, 3) = 2.d0*x5*d4_x5(1:3) + 2.d0*x6*d4_x6(1:3)
          d4_ads4(1:3, 4) = x1_2*d4_x3(1:3) + 2.d0*x2*x3*d4_x3(1:3)
          d4_ads4(1:3, 5) = (d4_x3(1:3)*x6 + x3*d4_x6(1:3))*(x1 + x2)
          d4_ads4(1:3, 6) = 2.d0*x2*x6*d4_x6(1:3) + 2.d0*x3*x5*d4_x5(1:3)
          d4_ads4(1:3, 7) = x2_2*d4_x6(1:3) + 2.d0*x3*x5*d4_x3(1:3) + x3_2*d4_x5(1:3)
          d4_ads4(1:3, 8) = x1*x4*d4_x5(1:3) + x2*x4*d4_x6(1:3) + d4_x3(1:3)*x5*x6 + d4_x5(1:3)*x3*x6 + d4_x6(1:3)*x5*x3
          d4_ads4(1:3, 9) = x4**2*d4_x6(1:3) + 2.d0*x4*x5*d4_x5(1:3) + d4_x5(1:3)*x6**2 + 2.d0*x5*x6*d4_x6(1:3)
          d4_ads4(1:3, 10) = d4_ads4(1:3, 2)*ads4(3) + ads4(2)*d4_ads4(1:3, 3)
          d4_ads4(1:3, 11) = 2.d0*ads4(2)*d4_ads4(1:3, 2)
          d4_ads4(1:3, 12) = 2.d0*ads4(3)*d4_ads4(1:3, 3)
        end if
        !else
        !   write (*,*) 'AAB, ABA, BAA NOT IMPLEMENTED'
        !   stop
        !   !*AAB
        !   if (i_type(i2)==i_type(i3)) then
        !      adp4(1) = x3
        !      adp4(2) = x1 + x3
        !      adp4(3) = x5 + x6
        !      adp4(4) = x4
        !      adp4(5) = x1**2 + x2**2
        !      adp4(6) = x5**2 + x6**2
        !
        !      s1 = 1.d0
        !      s2 = x1*x6 + x2*x5
        !   end if
        !
        !   !*BAA
        !   if (i_type(i3)==i_type(i4)) then
        !      adp4(1) = x1
        !      adp4(2) = x2 + x3
        !      adp4(3) = x4 + x5
        !      adp4(4) = x6
        !      adp4(5) = x2**2 + x3**2
        !      adp4(6) = x4**2 + x5**2
        !
        !      s1 = 1.d0
        !      s2 = x3*x4 + x2*x5
        !   end if
        !
        !   !*ABA
        !   if (i_type(i2)==i_type(i4)) then
        !      adp4(1) = x2
        !      adp4(2) = x1 + x3
        !      adp4(3) = x4 + x6
        !      adp4(4) = x5
        !      adp4(5) = x1**2 + x3**2
        !      adp4(6) = x4**2 + x6**2
        !
        !      s1 = 1.d0
        !      s2 = x1*x6 + x3*x4
        !   end if
        !
        !   if ( (i_type(i2)/=i_type(i3)).and.(i_type(i2)/=i_type(i4)).and.(i_type(i3)/=i_type(i4))) then
        !      adp4(1)= x1
        !      adp4(2)= x2
        !      adp4(3)= x3
        !      adp4(4)= x4
        !      adp4(5)= x5
        !      adp4(6)= x6
        !
        !      s1 = 1.d0
        !   end if
        !end if

        pmatrix(1:6, 0) = 1.d0
        d2_pmatrix(1:3, 1:6, 0) = 0.d0
        d3_pmatrix(1:3, 1:6, 0) = 0.d0
        d4_pmatrix(1:3, 1:6, 0) = 0.d0
        dd_pmatrix(1:6, 0) = 0.d0
        do im = 1, body_D_max(4)
          do ip = 1, 6
            pmatrix(ip, im) = adp4(ip)**im
            d2_pmatrix(1:3, ip, im) = dble(im)*adp4(ip)**(im - 1)*d2_adp4(1:3, ip)
            d3_pmatrix(1:3, ip, im) = dble(im)*adp4(ip)**(im - 1)*d3_adp4(1:3, ip)
            d4_pmatrix(1:3, ip, im) = dble(im)*adp4(ip)**(im - 1)*d4_adp4(1:3, ip)
            dd_pmatrix(ip, im) = dble(im)*adp4(ip)**(im - 1)
          end do
        end do

        tmp_cut = r_fcut(i2)*r_fcut(i3)*r_fcut(i4)
        if (desc_forces_bond) then
          d2_tmp_cut(1:3) = d_r_fcut(1:3, i2)*r_fcut(i3)*r_fcut(i4)
          d3_tmp_cut(1:3) = r_fcut(i2)*d_r_fcut(1:3, i3)*r_fcut(i4)
          d4_tmp_cut(1:3) = r_fcut(i2)*r_fcut(i3)*d_r_fcut(1:3, i4)
        end if

        if (debug_time) then
          time2 = MY_MPI_WTIME()
          time_invariants(4) = time_invariants(4) + time2 - time1
        end if

        do id = 1, dim_desc_body(4)
          iidk => kk4(1:6, id)

          pm1 = pmatrix(1, iidk(1))
          pm2 = pmatrix(2, iidk(2))
          pm3 = pmatrix(3, iidk(3))
          pm4 = pmatrix(4, iidk(4))
          pm5 = pmatrix(5, iidk(5))
          pm6 = pmatrix(6, iidk(6))
          ! tmp_pol = pmatrix(1,iidk(1))* pmatrix(2,iidk(2))*pmatrix(3,iidk(3))*pmatrix(4,iidk(4))*pmatrix(5,iidk(5))*pmatrix(6,iidk(6))
          pm123 = pm1*pm2*pm3
          pm456 = pm4*pm5*pm6
          ! tmp_pol = pm1 * pm2 * pm3 * pm4 * pm5 * pm6
          tmp_pol = pm123*pm456
          tmp_real(id) = tmp_real(id) + tmp_cut*ads4(bb4(id))*tmp_pol
          if (desc_forces_bond) then
            !p1=                     pmatrix(2,iidk(2))*pmatrix(3,iidk(3))*pmatrix(4,iidk(4))*pmatrix(5,iidk(5))*pmatrix(6,iidk(6))
            !p2= pmatrix(1,iidk(1))*                    pmatrix(3,iidk(3))*pmatrix(4,iidk(4))*pmatrix(5,iidk(5))*pmatrix(6,iidk(6))
            !p3= pmatrix(1,iidk(1))* pmatrix(2,iidk(2))*                   pmatrix(4,iidk(4))*pmatrix(5,iidk(5))*pmatrix(6,iidk(6))
            !p4= pmatrix(1,iidk(1))* pmatrix(2,iidk(2))*pmatrix(3,iidk(3))*                   pmatrix(5,iidk(5))*pmatrix(6,iidk(6))
            !p5= pmatrix(1,iidk(1))* pmatrix(2,iidk(2))*pmatrix(3,iidk(3))*pmatrix(4,iidk(4))*                   pmatrix(6,iidk(6))
            !p6= pmatrix(1,iidk(1))* pmatrix(2,iidk(2))*pmatrix(3,iidk(3))*pmatrix(4,iidk(4))*pmatrix(5,iidk(5))

            !p1 =       pm2 * pm3 * pm4 * pm5 * pm6
            !p2 = pm1 *       pm3 * pm4 * pm5 * pm6
            !p3 = pm1 * pm2 *       pm4 * pm5 * pm6
            !p4 = pm1 * pm2 * pm3 *       pm5 * pm6
            !p5 = pm1 * pm2 * pm3 * pm4 *       pm6
            !p6 = pm1 * pm2 * pm3 * pm4 * pm5

            p1 = pm2*pm3*pm456
            p2 = pm1*pm3*pm456
            p3 = pm1*pm2*pm456
            p4 = pm123*pm5*pm6
            p5 = pm123*pm4*pm6
            p6 = pm123*pm4*pm5


            d1 = dd_pmatrix(1, iidk(1))*p1
            d2 = dd_pmatrix(2, iidk(2))*p2
            d3 = dd_pmatrix(3, iidk(3))*p3
            d4 = dd_pmatrix(4, iidk(4))*p4
            d5 = dd_pmatrix(5, iidk(5))*p5
            d6 = dd_pmatrix(6, iidk(6))*p6
            !d_vect=(/  dd_pmatrix(1,iidk(1))*p1, &
            !           dd_pmatrix(2,iidk(2))*p2, &
            !           dd_pmatrix(3,iidk(3))*p3, &
            !           dd_pmatrix(4,iidk(4))*p4, &
            !           dd_pmatrix(5,iidk(5))*p5, &
            !           dd_pmatrix(6,iidk(6))*p6  /)
            !d2_tmp_pol(1:3) = d2_pmatrix(1:3,1,iidk(1))*p1 + &
            !                  d2_pmatrix(1:3,2,iidk(2))*p2 + &
            !                  d2_pmatrix(1:3,3,iidk(3))*p3 + &
            !                  d2_pmatrix(1:3,4,iidk(4))*p4 + &
            !                  d2_pmatrix(1:3,5,iidk(5))*p5 + &
            !                  d2_pmatrix(1:3,6,iidk(6))*p6
            !

            d2_tmp_pol(1:3) = d1*d2_adp4(1:3, 1) + &
                              d2*d2_adp4(1:3, 2) + &
                              d3*d2_adp4(1:3, 3) + &
                              d4*d2_adp4(1:3, 4) + &
                              d5*d2_adp4(1:3, 5) + &
                              d6*d2_adp4(1:3, 6)
            !do ix=1,3
            !d2_tmp_pol(ix) = ddot(6,d2_adp4(ix,:),1, d_vect(:),1)
            !end do
            !call dgemv('N',3,6,1.d0,d2_adp4,3,d_vect,1,0.d0,d2_tmp_pol,1)
            !call gemv(d2_adp4,d_vect,d2_tmp_pol)

            at = ads4(bb4(id))*tmp_pol
            tt = tmp_cut*tmp_pol
            ta = tmp_cut*ads4(bb4(id))
            !d2_tmp_real(id,i2,1:3) =  d2_tmp_real(id,i2,1:3) + d2_tmp_cut(1:3)*   ads4(bb4(id))    *   tmp_pol  + &
            !                                                      tmp_cut     *d2_ads4(1:3,bb4(id))*   tmp_pol  + &
            !                                                      tmp_cut     *   ads4(bb4(id))    *d2_tmp_pol(1:3)

            !d2_tmp_real(id,i2,1:3) =  d2_tmp_real(id,i2,1:3) + d2_tmp_cut(1:3)*   at  + &
            !                                                      tt     *d2_ads4(1:3,bb4(id)) + &
            !                                                      ta   *d2_tmp_pol(1:3)


            !d3_tmp_pol(1:3) = d3_pmatrix(1:3,1,iidk(1))*p1 + &
            !                  d3_pmatrix(1:3,2,iidk(2))*p2 + &
            !                  d3_pmatrix(1:3,3,iidk(3))*p3 + &
            !                  d3_pmatrix(1:3,4,iidk(4))*p4 + &
            !                  d3_pmatrix(1:3,5,iidk(5))*p5 + &
            !                  d3_pmatrix(1:3,6,iidk(6))*p6

            d3_tmp_pol(1:3) = d1*d3_adp4(1:3, 1) + &
                              d2*d3_adp4(1:3, 2) + &
                              d3*d3_adp4(1:3, 3) + &
                              d4*d3_adp4(1:3, 4) + &
                              d5*d3_adp4(1:3, 5) + &
                              d6*d3_adp4(1:3, 6)
            !do ix=1,3
            !d3_tmp_pol(ix) = ddot(6,d3_adp4(ix,:),1, d_vect(:),1)
            !end do
            !call dgemv('N',3,6,1.d0,d3_adp4,3,d_vect,1,0.d0,d3_tmp_pol,1)
            !call gemv(d3_adp4,d_vect,d3_tmp_pol)
            !d3_tmp_real(id,i3,1:3) =  d3_tmp_real(id,i3,1:3) + d3_tmp_cut(1:3)*   ads4(bb4(id))    *   tmp_pol  + &
            !                                                      tmp_cut     *d3_ads4(1:3,bb4(id))*   tmp_pol  + &
            !                                                      tmp_cut     *   ads4(bb4(id))    *d3_tmp_pol(1:3)
            !d3_tmp_real(id,i3,1:3) =  d3_tmp_real(id,i3,1:3) + d3_tmp_cut(1:3)*   at + &
            !                                                      tt    *d3_ads4(1:3,bb4(id))  + &
            !                                                      ta    *d3_tmp_pol(1:3)
            !d4_tmp_pol(1:3) = d4_pmatrix(1:3,1,iidk(1))*p1 + &
            !                  d4_pmatrix(1:3,2,iidk(2))*p2 + &
            !                  d4_pmatrix(1:3,3,iidk(3))*p3 + &
            !                  d4_pmatrix(1:3,4,iidk(4))*p4 + &
            !                  d4_pmatrix(1:3,5,iidk(5))*p5 + &
            !                  d4_pmatrix(1:3,6,iidk(6))*p6

            d4_tmp_pol(1:3) = d1*d4_adp4(1:3, 1) + &
                              d2*d4_adp4(1:3, 2) + &
                              d3*d4_adp4(1:3, 3) + &
                              d4*d4_adp4(1:3, 4) + &
                              d5*d4_adp4(1:3, 5) + &
                              d6*d4_adp4(1:3, 6)
            !do ix=1,3
            !d4_tmp_pol(ix) = ddot(6,d4_adp4(ix,:),1, d_vect(:),1)
            !end do
            !call dgemv('N',3,6,1.d0,d4_adp4,3,d_vect,1,0.d0,d4_tmp_pol,1)
            !call gemv(d4_adp4,d_vect,d4_tmp_pol)
            !d4_tmp_real(id,i4,1:3) =  d4_tmp_real(id,i4,1:3) + d4_tmp_cut(1:3)*   ads4(bb4(id))    *   tmp_pol  + &
            !                                                      tmp_cut     *d4_ads4(1:3,bb4(id))*   tmp_pol  + &
            !                                                      tmp_cut     *   ads4(bb4(id))    *d4_tmp_pol(1:3)
            do ix = 1, 3
              tmp_i2(id, ix) = d2_tmp_cut(ix)*at + tt*d2_ads4(ix, bb4(id)) + ta*d2_tmp_pol(ix)
              tmp_i3(id, ix) = d3_tmp_cut(ix)*at + tt*d3_ads4(ix, bb4(id)) + ta*d3_tmp_pol(ix)
              tmp_i4(id, ix) = d4_tmp_cut(ix)*at + tt*d4_ads4(ix, bb4(id)) + ta*d4_tmp_pol(ix)
            end do
            !d4_tmp_real(id,i4,1:3) =  d4_tmp_real(id,i4,1:3) + d4_tmp_cut(1:3)*   at  + &
            !                                                      tt     *d4_ads4(1:3,bb4(id))  + &
            !                                                      ta    *d4_tmp_pol(1:3)
          end if
          !9.1 tmp_pol = PRODUCT(adp4(:)**basis4(id)%ka(:))

          !tmp_pol=1.d0
          !do im=1,6
          !  ip=kk4(im,id)
          !  ip= iidk(im)
          !  tmp_pol = tmp_pol * pmatrix(im,ip)
          !end do

          !1.2 tmp_pol = pmatrix(1,basis4(id)%ka(1))* pmatrix(2,basis4(id)%ka(2))*pmatrix(3,basis4(id)%ka(3))*pmatrix(4,basis4(id)%ka(4))*pmatrix(5,basis4(id)%ka(5))*pmatrix(6,basis4(id)%ka(6))

          !0.7 tmp_pol = pmatrix(1,kk4(1,id))* pmatrix(2,kk4(2,id))*pmatrix(3,kk4(3,id))*pmatrix(4,kk4(4,id))*pmatrix(5,kk4(5,id))*pmatrix(6,kk4(6,id))
          !4.5 tmp_pol = PRODUCT(pmatrix(1:6,iidk(1:6)))!* pmatrix(2,basis4(id)%ka(2))*pmatrix(3,basis4(id)%ka(3))*pmatrix(4,basis4(id)%ka(4))*pmatrix(5,basis4(id)%ka(5))*pmatrix(6,basis4(id)%ka(6))
          !0.1 tmp_pol = tmp_cut*tmp_cut*tmp_cut*tmp_cut*tmp_cut*tmp_cut
          !tmp_real(id) =  tmp_real(id) + tmp_cut*ads4(basis4(id)%b)*tmp_pol

        end do                  ! dim_desc
        do ix = 1, 3
          call daxpy(dim_desc_body(4), 1.d0, tmp_i2(:, ix), 1, d2_tmp_real(:, i2, ix), 1)
          call daxpy(dim_desc_body(4), 1.d0, tmp_i3(:, ix), 1, d3_tmp_real(:, i3, ix), 1)
          call daxpy(dim_desc_body(4), 1.d0, tmp_i4(:, ix), 1, d4_tmp_real(:, i4, ix), 1)
          !   do ii=1, dim_desc_body(4)
          !      d2_tmp_real(ii, i2,ix) = d2_tmp_real(ii,i2,ix) + tmp_i2(ii,ix)
          !      d3_tmp_real(ii, i3,ix) = d3_tmp_real(ii,i3,ix) + tmp_i3(ii,ix)
          !      d4_tmp_real(ii, i4,ix) = d4_tmp_real(ii,i4,ix) + tmp_i4(ii,ix)
          !   end do
        end do

        if (debug_time) then
          time3 = MY_MPI_WTIME()
          time_polynomial(4) = time_polynomial(4) + time3 - time2
        end if

      end do                  ! i4
    end do                  ! i3
  end do                  ! i2

  if (desc_forces_bond) then
    tmp_force(1:dim_desc_body(4), 1:max_neigh_local, 1:3) = (d2_tmp_real(1:dim_desc_body(4), 1:max_neigh_local, 1:3) + &
                                                             d3_tmp_real(1:dim_desc_body(4), 1:max_neigh_local, 1:3) + &
                                                             d4_tmp_real(1:dim_desc_body(4), 1:max_neigh_local, 1:3))/4.d0
  end if
  tmp_real(:) = tmp_real(:)/4.d0                   ! because i2,i3,i4 are symmetric indexes
end subroutine bond_order_4



subroutine init_body()

  use ml_in_ndm_module, only: rangml
  use module_body_desc, only: deg_p2, deg_s2, &
                              deg_p3, deg_s3, deg_p4, deg_s4, &
                              dim_desc_body, body_D_max, l_body_order, &
                              basis, basis2, basis3, basis4, bb4, kk4
  implicit none

  integer  :: i, itot, itemp, ival
  integer  :: ka1, ka2, ka3, ka4, ka5, ka6, dim_k, m


  dim_desc_body(:) = 0
  ! INIT
  itot = 0

  ! order 2
  if (l_body_order(2)) then
    deg_p2(1) = 1
    deg_s2(1) = 0
    ! number of k for order n=2 : 2(2-1)/2 = 1
    dim_k = 1
    itemp = 0
    do i = 1, size(deg_s2)
      do m = 1, body_D_max(2)
        do ka1 = 0, body_D_max(2)
          ival = deg_s2(i) + ka1*deg_p2(1)
          if (ival /= m) cycle
          itemp = itemp + 1
          itot = itot + 1
        end do
      end do
    end do
    dim_desc_body(2) = itemp
  end if

  ! order 3
  if (l_body_order(3)) then
    deg_p3(1:3) = (/1, 2, 1/)
    deg_s3(1) = 0
    ! number of k for order n=3 : 3(3-1)/2 = 3
    dim_k = 3
    itemp = 0
    do i = 1, size(deg_s3)
      do m = 1, body_D_max(3)
        do ka1 = 0, body_D_max(3)
        do ka2 = 0, body_D_max(3)
        do ka3 = 0, body_D_max(3)
          ival = deg_s3(i) + ka1*deg_p3(1) + ka2*deg_p3(2) + ka3*deg_p3(3)
          if (ival /= m) cycle
          itemp = itemp + 1
          itot = itot + 1
        end do
        end do
        end do
      end do
    end do
    dim_desc_body(3) = itemp
  end if

  ! order 4
  if (l_body_order(4)) then
    deg_p4(1:6) = (/1, 1, 2, 2, 2, 3/)
    deg_s4(1:12) = (/0, 2, 2, 3, 3, 3, 3, 3, 3, 4, 4, 6/)
    itemp = 0
    ! number of k for order n=4 : 4(4-1)/2 = 6
    dim_k = 6
    do i = 1, size(deg_s4)
      do m = 1, body_D_max(4)
        do ka1 = 0, body_D_max(4)
        do ka2 = 0, body_D_max(4)
        do ka3 = 0, body_D_max(4)
        do ka4 = 0, body_D_max(4)
        do ka5 = 0, body_D_max(4)
        do ka6 = 0, body_D_max(4)
          ival = deg_s4(i) + ka1*deg_p4(1) + ka2*deg_p4(2) + ka3*deg_p4(3) &
                 + ka4*deg_p4(4) + ka5*deg_p4(5) + ka6*deg_p4(6)
          if (ival /= m) cycle
          itemp = itemp + 1
          itot = itot + 1
        end do
        end do
        end do
        end do
        end do
        end do
      end do
    end do
    dim_desc_body(4) = itemp
  end if


  if (l_body_order(5)) then
    if (rangml == 0) write (6, *) 'ML: FATAL ERROR init_body. Body order 5 not yet implemented'
    stop 'body order not implemented in init_boby'
  end if

  do i = 1, 5
    if (rangml == 0) then
      if (l_body_order(i)) write (6, '("ML: dimension from body order ", i4, " is .....................:",i7)') i, dim_desc_body(i)
    end if
  end do

  if (rangml == 0) then
    write (6, '("ML: Total dimension for this many-body expansion  is ......:",i7)') SUM(dim_desc_body(:))
    if (SUM(dim_desc_body(:)) == 0) then
      write (6, *) 'At least one many-body order should be activated in order to have the total dimension > 0'
      stop ' Total dimension for many-body descriptor is 0 '
    end if
  end if

  ! ALLOCATE
  if (allocated(basis)) deallocate (basis); allocate (basis(SUM(dim_desc_body(:))))

  itot = 0
  ! order 2
  if (l_body_order(2)) then
    if (allocated(basis2)) deallocate (basis2); allocate (basis2(dim_desc_body(2)))
    ! number of k for order n=2 : 2(2-1)/2 = 1
    dim_k = 1
    itemp = 0
    do i = 1, size(deg_s2)
      do m = 1, body_D_max(2)
        do ka1 = 0, body_D_max(2)
          ival = deg_s2(i) + ka1*deg_p2(1)
          if (ival /= m) cycle
          itemp = itemp + 1
          itot = itot + 1
          if (allocated(basis(itot)%ka)) deallocate (basis(itot)%ka); allocate (basis(itot)%ka(dim_k))
          basis(itot)%n = 2
          basis(itot)%b = i
          basis(itot)%Dmax = body_D_max(2)
          basis(itot)%Dcurrent = m
          basis(itot)%i_in_body = itemp
          basis(itot)%ka(1) = ka1
          if (allocated(basis2(itemp)%ka)) deallocate (basis2(itemp)%ka); allocate (basis2(itemp)%ka(dim_k))
          basis2(itemp)%n = 2
          basis2(itemp)%b = i
          basis2(itemp)%Dmax = body_D_max(2)
          basis2(itemp)%Dcurrent = m
          basis2(itemp)%i_in_body = itemp
          basis2(itemp)%ka(1) = ka1
        end do
      end do
    end do
  end if

  ! order 3
  if (l_body_order(3)) then
    if (allocated(basis3)) deallocate (basis3); allocate (basis3(dim_desc_body(3)))
    ! number of k for order n=3 : 3(3-1)/2 = 3
    dim_k = 3
    itemp = 0
    do i = 1, size(deg_s3)
      do m = 1, body_D_max(3)
        do ka1 = 0, body_D_max(3)
        do ka2 = 0, body_D_max(3)
        do ka3 = 0, body_D_max(3)
          ival = deg_s3(i) + ka1*deg_p3(1) + ka2*deg_p3(2) + ka3*deg_p3(3)
          if (ival /= m) cycle
          itemp = itemp + 1
          itot = itot + 1
          if (allocated(basis(itot)%ka)) deallocate (basis(itot)%ka); allocate (basis(itot)%ka(dim_k))
          basis(itot)%n = 3
          basis(itot)%b = i
          basis(itot)%Dmax = body_D_max(3)
          basis(itot)%Dcurrent = m
          basis(itot)%i_in_body = itemp
          basis(itot)%ka(1:3) = (/ka1, ka2, ka3/)
          if (allocated(basis3(itemp)%ka)) deallocate (basis3(itemp)%ka); allocate (basis3(itemp)%ka(dim_k))
          basis3(itemp)%n = 3
          basis3(itemp)%b = i
          basis3(itemp)%Dmax = body_D_max(3)
          basis3(itemp)%Dcurrent = m
          basis3(itemp)%i_in_body = itemp
          basis3(itemp)%ka(1:3) = (/ka1, ka2, ka3/)
        end do
        end do
        end do
      end do
    end do
  end if

  ! order 4
  if (l_body_order(4)) then
    if (allocated(basis4)) deallocate (basis4); allocate (basis4(dim_desc_body(4)))
    if (allocated(bb4)) deallocate (bb4); allocate (bb4(dim_desc_body(4)))
    if (allocated(kk4)) deallocate (kk4); allocate (kk4(6, dim_desc_body(4)))
    itemp = 0
    ! number of k for order n=4 : 4(4-1)/2 = 6
    dim_k = 6
    do i = 1, size(deg_s4)
      do m = 1, body_D_max(4)
        do ka1 = 0, body_D_max(4)
        do ka2 = 0, body_D_max(4)
        do ka3 = 0, body_D_max(4)
        do ka4 = 0, body_D_max(4)
        do ka5 = 0, body_D_max(4)
        do ka6 = 0, body_D_max(4)
          ival = deg_s4(i) + ka1*deg_p4(1) + ka2*deg_p4(2) + ka3*deg_p4(3) &
                 + ka4*deg_p4(4) + ka5*deg_p4(5) + ka6*deg_p4(6)
          if (ival /= m) cycle
          itemp = itemp + 1
          itot = itot + 1
          if (allocated(basis(itot)%ka)) deallocate (basis(itot)%ka); allocate (basis(itot)%ka(dim_k))
          basis(itot)%n = 4
          basis(itot)%b = i
          basis(itot)%Dmax = body_D_max(4)
          basis(itot)%Dcurrent = m
          basis(itot)%i_in_body = itemp
          basis(itot)%ka(1:6) = (/ka1, ka2, ka3, ka4, ka5, ka6/)
          bb4(itemp) = i
          if (allocated(basis4(itemp)%ka)) deallocate (basis4(itemp)%ka); allocate (basis4(itemp)%ka(dim_k))
          basis4(itemp)%n = 4
          basis4(itemp)%b = i
          basis4(itemp)%Dmax = body_D_max(4)
          basis4(itemp)%Dcurrent = m
          basis4(itemp)%i_in_body = itemp
          basis4(itemp)%ka(1:6) = (/ka1, ka2, ka3, ka4, ka5, ka6/)
          kk4(1:6, itemp) = (/ka1, ka2, ka3, ka4, ka5, ka6/)
        end do
        end do
        end do
        end do
        end do
        end do
      end do
    end do
  end if

end subroutine init_body
