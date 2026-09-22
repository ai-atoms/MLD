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

module compute_g2_mod
  USE module_kind_variables, ONLY: double
  implicit none

contains

  subroutine compute_g2(i_start_at, i_final_at, local_d_n_neigh, local_d_kind_neigh, iconf)

#ifdef MLD_NDM
    use gen_com_m, ONLY: A2cm, lperiod
    use gen_com_m_ml, ONLY: imm, bg, at
    use tab_imm_m_ml, ONLY: xp
#else
    use ondm_gen_com_m, ONLY: imm, A2cm, lperiod, bg, at, indi2
    use ondm_tab_imm_m, ONLY: iwmax2, xp
#endif
    use ml_in_ndm_module, ONLY: weighted, weighted_3ch, &
                                one_pi, imm_neigh, &
                                g2_eta, g2_rs, &
                                g2_dim, desc_forces, weighted, linvisible
    use derived_types, only: config_real, config_desc
    use module_neigh_local, only: r_cut 

#ifdef MLD_NDM
  use notperiod_mod, only: notperiod
#else
    use ondm_transform_coord, only: ondm_notperiod
#endif



    integer, intent(in)  :: i_start_at, i_final_at
    integer, dimension(imm), intent(out)   :: local_d_n_neigh
    integer, dimension(imm, imm_neigh), intent(out)    :: local_d_kind_neigh
    !double precision,dimension(g2_dim,imm),intent(out) :: local_g2_out
    !double precision,dimension(g2_dim,imm, 0:imm_neigh,3),intent(out) :: local_g2_deriv_out
    integer, optional, intent(in)    :: iconf

    real(double), dimension(:, :), allocatable   :: xpnp
    real(double), dimension(3) :: dxp_ji, ds
    integer  :: iw, iw1, iw2

    integer  :: p
    real(double)   :: fcut, r_ji, dfcut, fcut_w, dfcut_w, fcut_w_3ch, dfcut_w_3ch
    real(double)   :: factor_ia, factor_ja, factor_ia_3ch, factor_ja_3ch, g2_func

    integer  :: ia, ja, ia_n
    logical  :: small, desc_forces_local

    if ((i_start_at == 0) .and. (i_final_at == 0)) then
      local_d_kind_neigh(:, :) = 0
      local_d_n_neigh(:) = 0
      config_desc(iconf)%energy(:, :) = 0.d0
      !local_g2_deriv_out(:,:,:,:)=0.d0
      return
    end if
    desc_forces_local = desc_forces .and. (config_real(iconf)%has_force .or. config_real(iconf)%has_stress)
    small = .false.
    if (present(iconf)) then
      small = config_real(iconf)%small
    end if

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


    local_d_kind_neigh(:, :) = 0
    local_d_n_neigh(:) = 0
    config_desc(iconf)%energy(:, :) = 0d0
    if (desc_forces_local) config_desc(iconf)%force(:, :, :, :) = 0d0

#ifdef MLD_NDM
    iw2=0
!!$    if (i_start_at == 1) iw2 = 0
!!$    if (i_start_at > 1) iw2 = iwmax2(i_start_at - 1)
#else
    if (i_start_at == 1) iw2 = 0
    if (i_start_at > 1) iw2 = iwmax2(i_start_at - 1)
#endif

    do ja = i_start_at, i_final_at

      if (linvisible .and. weighted) then
        if (config_real(iconf)%invisible_per_type(config_real(iconf)%itype(ja))) cycle
      end if

      !begin small box or not 1/
#ifdef MLD_NDM
!!$      if (small) then
#else
      if (small) then
#endif
        iw1 = 1
        iw2 = config_real(iconf)%n_neigh(ja)
        !config_desc(iconf)%incell(ja,:)=.false.
#ifdef MLD_NDM
!!$      else
!!$        iw1 = iw2 + 1
!!$        iw2 = iwmax2(ja)
!!$      end if
#else
      else
        iw1 = iw2 + 1
        iw2 = iwmax2(ja)
      end if
#endif

      if (weighted) then
        factor_ja = config_real(iconf)%weight_per_type(config_real(iconf)%itype(ja))
        if (weighted_3ch) factor_ja_3ch = config_real(iconf)%weight_per_type_3ch(config_real(iconf)%itype(ja))
      else
        factor_ja = 1.d0
      end if

      ia_n = 0
      do iw = iw1, iw2
        !begin small box or not 2/
#ifdef MLD_NDM
!!$        if (small) then
#else
        if (small) then
#endif
          ia = config_real(iconf)%kind_neigh(ja, iw)
#ifdef MLD_NDM
!!$        else
!!$          ia = indi2(iw)
!!$          if (ja == ia) cycle
!!$        end if
#else
        else
          ia = indi2(iw)
          if (ja == ia) cycle
        end if
#endif


        if (linvisible .and. weighted) then
          if (config_real(iconf)%invisible_per_type(config_real(iconf)%itype(ia))) cycle
        end if

        !end   small box or not 2/
        if (weighted) then
          factor_ia = config_real(iconf)%weight_per_type(config_real(iconf)%itype(ia))
          if (weighted_3ch) factor_ia_3ch = config_real(iconf)%weight_per_type_3ch(config_real(iconf)%itype(ia))
        else
          factor_ia = 1.d0
        end if

#ifdef MLD_NDM
!!$        if (small) then
#else
        if (small) then
#endif
          r_ji = config_real(iconf)%r_ij(ja, iw)
          dxp_ji(:) = config_real(iconf)%u_ij(ja, iw, :)
          !incell = config_real(iconf)%incell(ja,iw)
#ifdef MLD_NDM
!!$        else
!!$          dxp_ji(1:3) = xpnp(1:3, ia) - xpnp(1:3, ja)
!!$          ds = MatMul(dxp_ji, bg)
!!$          WHERE ((ds .GT. 0.5d0) .OR. (ds .LT. -0.5d0))
!!$            ds(1:3) = ds(1:3) - Dble(Nint(ds(1:3)))
!!$          END WHERE
!!$          dxp_ji = MatMul(at, ds)/A2cm
!!$          r_ji = dsqrt(Sum(dxp_ji(1:3)**2))
!!$        end if
#else
        else
          dxp_ji(1:3) = xpnp(1:3, ia) - xpnp(1:3, ja)
          ds = MatMul(dxp_ji, bg)
          WHERE ((ds .GT. 0.5d0) .OR. (ds .LT. -0.5d0))
            ds(1:3) = ds(1:3) - Dble(Nint(ds(1:3)))
          END WHERE
          dxp_ji = MatMul(at, ds)/A2cm
          r_ji = dsqrt(Sum(dxp_ji(1:3)**2))
        end if
#endif

        if (r_ji >= r_cut) cycle
        ia_n = ia_n + 1

        fcut = 0.5d0*(dcos(one_pi*r_ji/r_cut) + 1.d0)
        if (desc_forces_local) dfcut = -0.5d0*one_pi*dsin(one_pi*r_ji/r_cut)/r_cut

        if (weighted) then
          fcut_w = fcut*factor_ia
          dfcut_w = dfcut*factor_ia
          if (weighted_3ch) then
            fcut_w_3ch = fcut*factor_ia_3ch
            dfcut_w_3ch = dfcut*factor_ia_3ch
          end if
        end if
        !debug write (*,*) fcut, dfcut, factor_ia
        do p = 1, g2_dim
          g2_func = dexp(-g2_eta(p)*(r_ji - g2_rs(p))**2)
          config_desc(iconf)%energy(p, ja) = config_desc(iconf)%energy(p, ja) + fcut*g2_func
          if (weighted) then
            config_desc(iconf)%energy(p + g2_dim, ja) = config_desc(iconf)%energy(p + g2_dim, ja) + fcut_w*g2_func
            if (weighted_3ch) config_desc(iconf)%energy(p + 2*g2_dim, ja) = config_desc(iconf)%energy(p + 2*g2_dim, ja) + fcut_w_3ch*g2_func
          end if
          if (desc_forces_local) then
            !if (.not.(allocated(config_desc(iconf)%force))) write (6,*) "FORCE NOT ALLOCATED", iconf, p , ja, ia_n
            config_desc(iconf)%force(p, ja, ia_n, 1:3) = (-2.d0*g2_eta(p)*(r_ji - g2_rs(p))*fcut + dfcut)*g2_func*dxp_ji(1:3)/r_ji
            if (weighted) then
              config_desc(iconf)%force(p + g2_dim, ja, ia_n, 1:3) = &
                (-2.d0*g2_eta(p)*(r_ji - g2_rs(p))*fcut_w + dfcut_w)* &
                g2_func*dxp_ji(1:3)/r_ji
              if (weighted_3ch) config_desc(iconf)%force(p + 2*g2_dim, ja, ia_n, 1:3) = &
                (-2.d0*g2_eta(p)*(r_ji - g2_rs(p))*fcut_w_3ch + dfcut_w_3ch)* &
                g2_func*dxp_ji(1:3)/r_ji
            end if
            config_desc(iconf)%force(p, ja, 0, 1:3) = config_desc(iconf)%force(p, ja, 0, 1:3) - config_desc(iconf)%force(p, ja, ia_n, 1:3)
            if (weighted) then
              config_desc(iconf)%force(p + g2_dim, ja, 0, 1:3) = config_desc(iconf)%force(p + g2_dim, ja, 0, 1:3) - config_desc(iconf)%force(p + g2_dim, ja, ia_n, 1:3)
              if (weighted_3ch) config_desc(iconf)%force(p + 2*g2_dim, ja, 0, 1:3) = config_desc(iconf)%force(p + 2*g2_dim, ja, 0, 1:3) - config_desc(iconf)%force(p + g2_dim, ja, ia_n, 1:3)
            end if
          end if
        end do
        local_d_kind_neigh(ja, ia_n) = ia
      end do                  ! ia
      if (weighted) then
        config_desc(iconf)%energy(g2_dim + 1:2*g2_dim, :) = config_desc(iconf)%energy(g2_dim + 1:2*g2_dim, :)*factor_ja             ! /dble(ia_n)
        if (desc_forces_local) config_desc(iconf)%force(g2_dim + 1:2*g2_dim, :, :, :) = config_desc(iconf)%force(g2_dim + 1:2*g2_dim, :, :, :)*factor_ja     ! /dble(ia_n)
        if (weighted_3ch) then
          config_desc(iconf)%energy(2*g2_dim + 1:3*g2_dim, :) = config_desc(iconf)%energy(2*g2_dim + 1:3*g2_dim, :)*factor_ja_3ch     ! /dble(ia_n)
          if (desc_forces_local) config_desc(iconf)%force(2*g2_dim + 1:3*g2_dim, :, :, :) = config_desc(iconf)%force(2*g2_dim + 1:2*g2_dim, :, :, :)*factor_ja_3ch                      ! /dble(ia_n)
        end if

      end if
      local_d_n_neigh(ja) = ia_n
    end do                  ! ja

    deallocate (xpnp)

  end subroutine compute_g2


  subroutine init_g2()
    use ml_in_ndm_module, only: rangml, g2_dim, n_g2_eta, n_g2_rs, strict_behler, &
                                eta_max_g2, eta_min_g2, rs_min_g2, rs_max_g2, g2_eta, g2_rs
    use module_neigh_local, only: r_cut 
    use mld_logger
    real(kind(0.d0)), dimension(:), allocatable  :: vd1
    real(kind(0.d0)), parameter      :: b2a = 0.52917721067d0
    integer  :: icount, p_e, p_r

    g2_dim = (n_g2_eta)*n_g2_rs

    if (strict_behler) then
      g2_dim = 8
    end if

    if (allocated(g2_eta)) deallocate (g2_eta); allocate (g2_eta(g2_dim))
    if (allocated(g2_rs)) deallocate (g2_rs); allocate (g2_rs(g2_dim))

    if (strict_behler) then

      r_cut = 11.338d0*b2a
      if (allocated(vd1)) deallocate (vd1); allocate (vd1(g2_dim))
      vd1(:) = (/0.001d0, 0.010d0, 0.020d0, 0.035d0, 0.060d0, 0.100d0, 0.200d0, 0.400d0/)
      g2_eta(1:g2_dim) = vd1(1:g2_dim)/b2a**2
      g2_rs(:) = 0.d0

    else

      icount = 0
      do p_e = 1, n_g2_eta
        do p_r = 1, n_g2_rs
          icount = icount + 1
          !old_G2: g2_eta(icount) = 1.d-2 + dble(p_e - 1)*(eta_max_g2 - 1.d-2)/dble(n_g2_eta - 1)
          if (n_g2_eta==1) then
            g2_eta(icount) = eta_min_g2
          else  
            g2_eta(icount) = 1.d0/(eta_min_g2 + dble(p_e-1)*(eta_max_g2 - eta_min_g2)/dble(n_g2_eta-1))**2
          end if 
          !old G2: g2_rs(icount) = 0.d0 + dble(p_r - 1)
          if (n_g2_rs == 1 ) then 
            g2_rs(icount) = rs_min_g2
          else
            g2_rs(icount) = rs_min_g2 + dble(p_r-1)*(rs_max_g2-rs_min_g2)/dble(n_g2_rs-1)
          end if 
        end do
      end do

      if (icount /= g2_dim) then
        if (rangml == 0) write (6, *) 'Error in subroutine init_g2'
        stop 'dimension of g2 is incorrect'
      end if

    end if

  end subroutine init_g2

end module
