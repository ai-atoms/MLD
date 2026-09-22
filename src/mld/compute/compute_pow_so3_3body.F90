! HND XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
! HND X
! HND X   MiLaDy (Machine Learning Dynamics) 
! HND X   
! HND X
! HND X   Milady was written and designed by Mihai-Cosmin Marinica, Alexandra M. Goryaeva 
! HND X   Copyright 2015-2023.
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

module module_compute_pow_so3_3body
contains
subroutine compute_pow_so3_3body(i_start_at, i_final_at, d_n_neigh, d_kind_neigh, iconf)

#if(PARA)
  !use mpi
  use mld_mpi 
#endif
  USE module_kind_variables, ONLY: double
#ifdef MLD_NDM
  use gen_com_m, ONLY: A2cm, lperiod
  use gen_com_m_ml, ONLY: imm, bg, at
  use tab_imm_m_ml, ONLY: xp
#else
  use ondm_gen_com_m, ONLY: imm,  lperiod
  use ondm_tab_imm_m, ONLY: iwmax2, xp
#endif
  use derived_types, only: config_real, config_desc
  use ml_in_ndm_module, ONLY: rangml, imm_neigh, l_max, &
                              desc_forces, linvisible, weighted, weighted_3ch, &
                              one_pi
#ifdef MLD_NDM
  use notperiod_mod
#endif

  use module_neigh_local, only: r_central, i_central, i_type, i_type_db, tmp_dxp, tmp_xp, &
                                max_neigh_local, iw2, r_cut, build_local_neighbours_ja
  use time_check_general, only: time, tot_time, debug_time, MY_MPI_WTIME
  use module_so3, only: n_rbf_so3, pow_so3_dim, clmn, dclmn, ini_rbf_so3, &
                        clmn_w, dclmn_w, clmn_w_3ch, dclmn_w_3ch
  use compute_pow_so3_mod, only: spherical_3d
#ifndef MLD_NDM
  use ondm_transform_coord, only: ondm_notperiod
#endif
  implicit none

  integer, intent(in)  :: i_start_at, i_final_at
  integer, dimension(imm), intent(out)   :: d_n_neigh
  integer, dimension(imm, imm_neigh), intent(out)    :: d_kind_neigh
  integer, optional    :: iconf

  logical  :: small
  real(double), dimension(:, :), allocatable   :: xpnp
  real(double), dimension(3) :: dxp_ji
  integer  :: ia, ja, ia_n, i_count, n1, n2, i1, i2

  integer  :: l, m
  real(double)   :: r_ji, tmp_fact

  real(double), dimension(pow_so3_dim)   :: tmp_pow_so3_out
  real(double), dimension(pow_so3_dim, 0:imm_neigh, 3)     :: tmp_pow_so3_deriv_out

  real(double), dimension(pow_so3_dim)   :: tmp_pow_so3_out_w
  real(double), dimension(pow_so3_dim, 0:imm_neigh, 3)     :: tmp_pow_so3_deriv_out_w

  real(double), dimension(pow_so3_dim)   :: tmp_pow_so3_out_w_3ch
  real(double), dimension(pow_so3_dim, 0:imm_neigh, 3)     :: tmp_pow_so3_deriv_out_w_3ch

  real(double)   :: factor_ia, factor_ja, factor_ia_3ch, factor_ja_3ch
  logical  :: desc_forces_local

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

  d_n_neigh(:) = 0
  d_kind_neigh(:, :) = 0
  config_desc(iconf)%energy(:, :) = 0.d0
  if (desc_forces_local) config_desc(iconf)%force(:, :, :, :) = 0.d0
#ifdef MLD_NDM
iw2=0
!!$    if (i_start_at == 1) iw2 = 0
!!$    if (i_start_at > 1) iw2 = iwmax2(i_start_at - 1)
!!$
#else
  if (i_start_at == 1) iw2 = 0
  if (i_start_at > 1) iw2 = iwmax2(i_start_at - 1)
#endif

  if (debug_time) time(1) = MY_MPI_WTIME()
  do ja = i_start_at, i_final_at
    ! begin small box or not 1/

    if (debug_time) time(3) = MY_MPI_WTIME()
    if (linvisible .and. weighted) then
      if (config_real(iconf)%invisible_per_type(config_real(iconf)%itype(ja))) cycle
    end if

    factor_ja = 1.d0
    factor_ja_3ch = 1.d0
    if (weighted) then
      factor_ja = config_real(iconf)%weight_per_type(config_real(iconf)%itype(ja))
      if (weighted_3ch) factor_ja_3ch = config_real(iconf)%weight_per_type_3ch(config_real(iconf)%itype(ja))
    end if

    tmp_pow_so3_out(:) = 0.d0
    tmp_pow_so3_deriv_out(:, :, :) = 0.d0
    clmn(:, :) = (0.d0, 0.d0)
    if (weighted) then
      tmp_pow_so3_out_w(:) = 0.d0
      tmp_pow_so3_deriv_out_w(:, :, :) = 0.d0
      clmn_w(:, :) = (0.d0, 0.d0)
      if (weighted_3ch) then
        tmp_pow_so3_out_w_3ch(:) = 0.d0
        tmp_pow_so3_deriv_out_w_3ch(:, :, :) = 0.d0
        clmn_w_3ch(:, :) = (0.d0, 0.d0)
      end if
    end if
    ia_n = 0

    call build_local_neighbours_ja( iconf, ja, imm, xpnp, r_cut, d_n_neigh, d_kind_neigh, &
                                    r_central, i_type, i_type_db, i_central, tmp_dxp, tmp_xp, &
                                    max_neigh_local, iw2)

    if (debug_time) then
      time(4) = MY_MPI_WTIME()
      tot_time(2) = tot_time(2) + time(4) - time(3)
    end if

    do ia_n = 1, max_neigh_local
      ! begin small box or not 2/
      ia = i_central(ia_n)
      r_ji = r_central(ia_n)
      dxp_ji = tmp_dxp(:, ia_n)
      if (debug_time) time(7) = MY_MPI_WTIME()

      factor_ia = 1.d0
      factor_ia_3ch = 1.d0
      if (weighted) then
        if (linvisible) then
          if (config_real(iconf)%invisible_per_type(config_real(iconf)%itype(ia))) cycle
        end if
        factor_ia = config_real(iconf)%weight_per_type(config_real(iconf)%itype(ia))
        if (weighted_3ch)  factor_ia_3ch = config_real(iconf)%weight_per_type_3ch(config_real(iconf)%itype(ia))
      end if

      if (debug_time) then
        time(8) = MY_MPI_WTIME()
        tot_time(5) = tot_time(5) + time(8) - time(7)
      end if

      ! this gives clmn and d_clmn ...
      call spherical_3d(ia_n, factor_ia, factor_ia_3ch, desc_forces_local, dxp_ji, r_ji)

      if (debug_time) then
        time(9) = MY_MPI_WTIME()
        tot_time(6) = tot_time(6) + time(9) - time(8)
      end if
    end do                  ! ia end of neighbours iterations


    ! update the componenet from the central atom ....
    ! is working worst is many cases ... but should be correct.
    !  ----             DO NOT ERASE          ----   !
    i_count = 0
    do n1 = ini_rbf_so3, n_rbf_so3
      do l = 0, l_max
        i_count = (n1 - ini_rbf_so3)*(l_max + 1) + l + 1
        tmp_fact =  dsqrt( (2.d0*dble(l)+1) / (4.d0 * one_pi) )
        clmn(0, i_count) = clmn(0, i_count) + tmp_fact
        if (weighted) then
          clmn_w(0, i_count) = clmn_w(0, i_count) + tmp_fact*factor_ja
          if (weighted_3ch) clmn_w_3ch(0, i_count) = clmn_w_3ch(0, i_count) + tmp_fact*factor_ja_3ch
        end if
      end do
    end do
    !  ----             DO NOT ERASE          ----   !

    d_n_neigh(ja) = max_neigh_local
    if (debug_time) then
      time(5) = MY_MPI_WTIME()
      tot_time(3) = tot_time(3) + time(5) - time(4)
    end if
    i_count = 0
    do n1 = ini_rbf_so3, n_rbf_so3
      do n2 = ini_rbf_so3, n_rbf_so3
        do l = 0, l_max
          i_count = i_count + 1
          i1 = (n1 - ini_rbf_so3)*(l_max + 1) + l + 1
          i2 = (n2 - ini_rbf_so3)*(l_max + 1) + l + 1
          do m = -l, l
            tmp_pow_so3_out(i_count) = tmp_pow_so3_out(i_count) + real(dconjg(clmn(m, i1))*clmn(m, i2), kind(0.d0))
            do ia = 1, d_n_neigh(ja)
              if (desc_forces_local) tmp_pow_so3_deriv_out(i_count, ia, 1:3) = tmp_pow_so3_deriv_out(i_count, ia, 1:3) + &
                                                                               real(dconjg(clmn(m, i1))*dclmn(m, i2, ia, 1:3) + clmn(m, i2)*dconjg(dclmn(m, i1, ia, 1:3)), kind(0.d0))
            end do
            if (weighted) then
              tmp_pow_so3_out_w(i_count) = tmp_pow_so3_out_w(i_count) + real(dconjg(clmn_w(m, i1))*clmn_w(m, i2), kind(0.d0))
              do ia = 1, d_n_neigh(ja)
                if (desc_forces_local) tmp_pow_so3_deriv_out_w(i_count, ia, 1:3) = tmp_pow_so3_deriv_out_w(i_count, ia, 1:3) + &
                                                                               real(dconjg(clmn_w(m, i1))*dclmn_w(m, i2, ia, 1:3) + clmn_w(m, i2)*dconjg(dclmn_w(m, i1, ia, 1:3)), kind(0.d0))
              end do
              if (weighted_3ch) then
                tmp_pow_so3_out_w_3ch(i_count) = tmp_pow_so3_out_w_3ch(i_count) + real(dconjg(clmn_w_3ch(m, i1))*clmn_w_3ch(m, i2), kind(0.d0))
                do ia = 1, d_n_neigh(ja)
                  if (desc_forces_local) tmp_pow_so3_deriv_out_w_3ch(i_count, ia, 1:3) = tmp_pow_so3_deriv_out_w_3ch(i_count, ia, 1:3) + &
                                                                               real(dconjg(clmn_w_3ch(m, i1))*dclmn_w_3ch(m, i2, ia, 1:3) + clmn_w_3ch(m, i2)*dconjg(dclmn_w_3ch(m, i1, ia, 1:3)), kind(0.d0))
                end do
              end if
            end if

          end do
        end do                  ! m
      end do                  ! p
    end do                 ! l

    ! forces .......
    if (desc_forces_local) then
      do ia = 1, d_n_neigh(ja)
        config_desc(iconf)%force(1:pow_so3_dim, ja, ia, 1:3) = tmp_pow_so3_deriv_out(1:pow_so3_dim, ia, 1:3)
        config_desc(iconf)%force(1:pow_so3_dim, ja, 0, 1:3) = config_desc(iconf)%force(1:pow_so3_dim, ja, 0, 1:3) - tmp_pow_so3_deriv_out(1:pow_so3_dim, ia, 1:3)
      end do
      if (weighted) then
        do ia = 1, d_n_neigh(ja)
          config_desc(iconf)%force(pow_so3_dim+1:2*pow_so3_dim, ja, ia, 1:3) = tmp_pow_so3_deriv_out_w(1:pow_so3_dim, ia, 1:3)!*factor_ja
          config_desc(iconf)%force(pow_so3_dim+1:2*pow_so3_dim, ja, 0, 1:3) = config_desc(iconf)%force(pow_so3_dim+1:2*pow_so3_dim, ja, 0, 1:3) - tmp_pow_so3_deriv_out_w(1:pow_so3_dim, ia, 1:3)!*factor_ja
        end do
        if (weighted_3ch) then
          do ia = 1, d_n_neigh(ja)
            config_desc(iconf)%force(2*pow_so3_dim+1:3*pow_so3_dim, ja, ia, 1:3) = tmp_pow_so3_deriv_out_w_3ch(1:pow_so3_dim, ia, 1:3)!*factor_ja_3ch
            config_desc(iconf)%force(2*pow_so3_dim+1:3*pow_so3_dim, ja, 0, 1:3) = config_desc(iconf)%force(2*pow_so3_dim+1:3*pow_so3_dim, ja, 0, 1:3) - tmp_pow_so3_deriv_out_w_3ch(1:pow_so3_dim, ia, 1:3)!*factor_ja_3ch
          end do
        end if
      end if
    end if

    ! energy .......
    config_desc(iconf)%energy(1:pow_so3_dim, ja) = tmp_pow_so3_out(1:pow_so3_dim)
    if (weighted) then
      config_desc(iconf)%energy(pow_so3_dim+1:2*pow_so3_dim, ja) = tmp_pow_so3_out_w(1:pow_so3_dim)!*factor_ja
      if (weighted_3ch) then
        config_desc(iconf)%energy(2*pow_so3_dim+1:3*pow_so3_dim, ja) = tmp_pow_so3_out_w_3ch(1:pow_so3_dim)!*factor_ja_3ch
      end if
    end if

    if (debug_time) then
      time(6) = MY_MPI_WTIME()
      tot_time(4) = tot_time(4) + time(6) - time(5)
    end if

  end do                  ! ja

  if (debug_time) then
    time(2) = MY_MPI_WTIME()
    tot_time(1) = tot_time(1) + time(2) - time(1)
  end if


  if (debug_time) then
    call repport_time(2, 0.d0, tot_time(1), "ML: full")
    call repport_time(2, 0.d0, tot_time(2), "ML: main loop: init before nn  ")
    call repport_time(2, 0.d0, tot_time(3), "ML: main loop: nn loop ")
    call repport_time(2, 0.d0, tot_time(4), "ML: main loop: final ")
    call repport_time(4, 0.d0, tot_time(5), "ML: nn loop: nn search")
    call repport_time(4, 0.d0, tot_time(6), "ML: nn loop: cmm/dcmm")
    call repport_time(6, 0.d0, tot_time(8), "ML: spher 3D : radial")
    call repport_time(6, 0.d0, tot_time(9), "ML: spher 3D :  spher")
    call repport_time(6, 0.d0, tot_time(10), "ML: spher 3D :    cmm")
  end if

  deallocate (xpnp)

end subroutine compute_pow_so3_3body
end module module_compute_pow_so3_3body


subroutine init_pow_so3_3body
  use ml_in_ndm_module, only: l_max,  imm_neigh, &
                              weighted, weighted_3ch
  use module_so3, only: n_rbf_so3, pow_so3_dim, radial_pow_so3, radial_bartok, radial_sgg, &
                        clmn, dclmn, ini_rbf_so3, clmn_w, dclmn_w, &
                        clmn_w_3ch, dclmn_w_3ch

  implicit none

  if (radial_pow_so3 == radial_bartok) then
    ini_rbf_so3 = 1
    pow_so3_dim = int((1 + l_max))*n_rbf_so3**2
  end if


  if (radial_pow_so3 == radial_sgg) then
    ini_rbf_so3 = 0
    pow_so3_dim = int((1 + l_max))*(n_rbf_so3 + 1)**2
  end if

  if (allocated(clmn)) deallocate (clmn); allocate (clmn(-l_max:l_max, pow_so3_dim))
  if (allocated(dclmn)) deallocate (dclmn); allocate (dclmn(-l_max:l_max, pow_so3_dim, imm_neigh, 3))
  if (weighted) then
    if (allocated(clmn_w)) deallocate (clmn_w); allocate (clmn_w(-l_max:l_max, pow_so3_dim))
    if (allocated(dclmn_w)) deallocate (dclmn_w); allocate (dclmn_w(-l_max:l_max, pow_so3_dim, imm_neigh, 3))
    if (weighted_3ch) then
     if (allocated(clmn_w_3ch)) deallocate (clmn_w_3ch); allocate (clmn_w_3ch(-l_max:l_max, pow_so3_dim))
     if (allocated(dclmn_w_3ch)) deallocate (dclmn_w_3ch); allocate (dclmn_w_3ch(-l_max:l_max, pow_so3_dim, imm_neigh, 3))
    end if
  end if

end subroutine init_pow_so3_3body
