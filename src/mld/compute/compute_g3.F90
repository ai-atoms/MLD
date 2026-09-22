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

module compute_g3_mod
  USE module_kind_variables, ONLY: double

  implicit none

contains

   subroutine compute_g3(i_start_at, i_final_at, local_d_n_neigh, local_d_kind_neigh, iconf)

    USE module_kind_variables, ONLY: double
#ifdef MLD_NDM
    use gen_com_m, ONLY: A2cm, lperiod
    use gen_com_m_ml, ONLY: imm, bg, at
    use tab_imm_m_ml, ONLY: xp
#else
    use ondm_gen_com_m, ONLY: imm, A2cm, lperiod, bg, at, indi2
    use ondm_tab_imm_m, ONLY: iwmax2, xp
#endif
    use ml_in_ndm_module, ONLY: desc_forces, debug, rangml, weighted, &
                                one_pi, g3_eta, g3_zeta, g3_lambda, g3_dim, imm_neigh, linvisible

    use derived_types, only: config_real, config_desc
    use module_neigh_local, only: r_cut 

#ifdef MLD_NDM
  use notperiod_mod
#else
    use ondm_transform_coord, only: ondm_notperiod
#endif


    integer, intent(in)  :: i_start_at, i_final_at
    integer, dimension(imm), intent(out)   :: local_d_n_neigh
    integer, dimension(imm, imm_neigh), intent(out)    :: local_d_kind_neigh
    integer, optional, intent(in)    :: iconf

    real(double), dimension(:, :), allocatable   :: xpnp
    real(double), dimension(3) :: dxp_ji, dxp_jk, dxp_ik, ds
    integer  :: ia, ia_n, ja, ka, ka_n, iw, iw1, iw2

    integer  :: iz, p
    double precision     :: fcut_ji, fcut_jk, fcut_ik
    double precision     :: r_ji, r_jk, r_ik
    double precision     :: r2_ji, r2_jk, r2_ik

    double precision     :: cos_jik, ang, rad, dang
    double precision, dimension(3)   :: dfcut_ji, dfcut_jk, dfcut_ik, cdxp_ik, cdxp_ji, cdxp_jk
    double precision, dimension(3)   :: dcos_jik_ji, dcos_jik_jk
    double precision, dimension(3)   :: drad_ji, drad_jk, drad_ik
    double precision     :: factor_ia, factor_ja, factor_ka

    logical  :: small, desc_forces_local


    if ((i_start_at == 0) .and. (i_final_at == 0)) then
      local_d_kind_neigh(:, :) = 0
      local_d_n_neigh(:) = 0
      config_desc(iconf)%energy(:, :) = 0.d0
      !config_desc(iconf)%force(:,:,:,:)=0.d0
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
    !call cryst_to_cart (imm, xpnp, bg, -1)

    local_d_kind_neigh(:, :) = 0
    local_d_n_neigh(:) = 0
    config_desc(iconf)%energy(:, :) = 0.d0
    if (desc_forces_local) config_desc(iconf)%force(:, :, :, :) = 0.d0

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

      if (debug) then
        if (mod(ja - 1, 10) == 0) then
          if (rangml == 0) write (6, '("in bso4 i_start_at i_final_at  ja:  ",3i9)') i_start_at, i_final_at, ja
        end if
      end if
      if (weighted) then
        factor_ja = config_real(iconf)%weight_per_type(config_real(iconf)%itype(ja))
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
        !end   small box or not 2/

        if (linvisible .and. weighted) then
          if (config_real(iconf)%invisible_per_type(config_real(iconf)%itype(ia))) cycle
        end if

        if (weighted) then
          factor_ia = config_real(iconf)%weight_per_type(config_real(iconf)%itype(ia))
        else
          factor_ia = 1.d0
        end if

#ifdef MLD_NDM
!!$        if (small) then
#else
        if (small) then
#endif
          r_ji = config_real(iconf)%r_ij(ja, iw)
          r2_ji = r_ji**2
          dxp_ji(:) = config_real(iconf)%u_ij(ja, iw, :)
#ifdef MLD_NDM
!!$        else
!!$          dxp_ji(1:3) = xpnp(1:3, ia) - xpnp(1:3, ja)
!!$          ds = MatMul(dxp_ji, bg)
!!$          WHERE ((ds .GT. 0.5d0) .OR. (ds .LT. -0.5d0))
!!$            ds(1:3) = ds(1:3) - Dble(Nint(ds(1:3)))
!!$          END WHERE
!!$          dxp_ji = MatMul(at, ds)/A2cm
!!$          r2_ji = Sum(dxp_ji(1:3)**2)
!!$          r_ji = dsqrt(r2_ji)
!!$        end if
#else
        else
          dxp_ji(1:3) = xpnp(1:3, ia) - xpnp(1:3, ja)
          ds = MatMul(dxp_ji, bg)
          WHERE ((ds .GT. 0.5d0) .OR. (ds .LT. -0.5d0))
            ds(1:3) = ds(1:3) - Dble(Nint(ds(1:3)))
          END WHERE
          dxp_ji = MatMul(at, ds)/A2cm
          r2_ji = Sum(dxp_ji(1:3)**2)
          r_ji = dsqrt(r2_ji)
        end if
#endif
        if (r_ji >= r_cut) cycle
        ia_n = ia_n + 1
        !
        cdxp_ji(1:3) = dxp_ji(1:3)/r_ji
        fcut_ji = 0.5d0*(cos(one_pi*r_ji/r_cut) + 1.d0)*factor_ia*factor_ja
        if (desc_forces_local) dfcut_ji(1:3) = -0.5d0*one_pi*sin(one_pi*r_ji/r_cut)*cdxp_ji(1:3)/r_cut*factor_ia*factor_ja
        ka_n = 0
        do iz = iw1, iw2
          !begin small box or not 2/
#ifdef MLD_NDM
!!$          if (small) then
#else
          if (small) then
#endif
            ka = config_real(iconf)%kind_neigh(ja, iz)
#ifdef MLD_NDM
!!$          else
!!$            ka = indi2(iz)
!!$            if ((ja == ka) .or. (ia == ka)) cycle
!!$          end if
#else
          else
            ka = indi2(iz)
            if ((ja == ka) .or. (ia == ka)) cycle
          end if
#endif
          !end   small box or not 2/


          if (linvisible .and. weighted) then
            if (config_real(iconf)%invisible_per_type(config_real(iconf)%itype(ka))) cycle
          end if

          if (weighted) then
            factor_ka = config_real(iconf)%weight_per_type(config_real(iconf)%itype(ka))
          else
            factor_ka = 1.d0
          end if

#ifdef MLD_NDM
!!$          if (small) then
#else
          if (small) then
#endif
            r_jk = config_real(iconf)%r_ij(ja, iz)
            r2_jk = r_jk**2
            dxp_jk(:) = config_real(iconf)%u_ij(ja, iz, :)
            dxp_ik(:) = dxp_jk(:) - dxp_ji(:)
            r2_ik = Sum(dxp_ik(:)**2)
            r_ik = dsqrt(r2_ik)
#ifdef MLD_NDM
!!$          else
!!$            dxp_jk(1:3) = xpnp(1:3, ka) - xpnp(1:3, ja)
!!$            ds = MatMul(dxp_jk, bg)
!!$            WHERE ((ds .GT. 0.5d0) .OR. (ds .LT. -0.5d0))
!!$              ds(1:3) = ds(1:3) - Dble(Nint(ds(1:3)))
!!$            END WHERE
!!$            dxp_jk = MatMul(at, ds)/A2cm
!!$            r2_jk = Sum(dxp_jk(1:3)**2)
!!$            r_jk = dsqrt(r2_jk)
!!$
!!$            dxp_ik(1:3) = xpnp(1:3, ka) - xpnp(1:3, ia)
!!$            ds = MatMul(dxp_ik, bg)
!!$            WHERE ((ds .GT. 0.5d0) .OR. (ds .LT. -0.5d0))
!!$              ds(1:3) = ds(1:3) - Dble(Nint(ds(1:3)))
!!$            END WHERE
!!$            dxp_ik = MatMul(at, ds)/A2cm
!!$            r2_ik = Sum(dxp_ik(1:3)**2)
!!$            r_ik = dsqrt(r2_ik)
!!$          end if
#else
          else
            dxp_jk(1:3) = xpnp(1:3, ka) - xpnp(1:3, ja)
            ds = MatMul(dxp_jk, bg)
            WHERE ((ds .GT. 0.5d0) .OR. (ds .LT. -0.5d0))
              ds(1:3) = ds(1:3) - Dble(Nint(ds(1:3)))
            END WHERE
            dxp_jk = MatMul(at, ds)/A2cm
            r2_jk = Sum(dxp_jk(1:3)**2)
            r_jk = dsqrt(r2_jk)

            dxp_ik(1:3) = xpnp(1:3, ka) - xpnp(1:3, ia)
            ds = MatMul(dxp_ik, bg)
            WHERE ((ds .GT. 0.5d0) .OR. (ds .LT. -0.5d0))
              ds(1:3) = ds(1:3) - Dble(Nint(ds(1:3)))
            END WHERE
            dxp_ik = MatMul(at, ds)/A2cm
            r2_ik = Sum(dxp_ik(1:3)**2)
            r_ik = dsqrt(r2_ik)
          end if
#endif

          if (r_jk >= r_cut) cycle
          ka_n = ka_n + 1
          if (r_ik >= r_cut) cycle
          if (r_ik <= 1d-20) cycle


          fcut_jk = 0.5d0*(cos(one_pi*r_jk/r_cut) + 1d0)*factor_ja*factor_ka
          fcut_ik = 0.5d0*(cos(one_pi*r_ik/r_cut) + 1d0)*factor_ia*factor_ka

          if (desc_forces_local) then
            cdxp_ik(1:3) = dxp_ik(1:3)/r_ik
            cdxp_jk(1:3) = dxp_jk(1:3)/r_jk
            dfcut_jk(1:3) = -0.5d0*one_pi*sin(one_pi*r_jk/r_cut)*cdxp_jk(1:3)/r_cut*factor_ja*factor_ka
            dfcut_ik(1:3) = -0.5d0*one_pi*sin(one_pi*r_ik/r_cut)*cdxp_ik(1:3)/r_cut*factor_ia*factor_ka
          end if

          cos_jik = dot_product(dxp_ji(1:3)/r_ji, dxp_jk(1:3)/r_jk)

          if (desc_forces_local) then
            dcos_jik_ji(1:3) = (cdxp_jk(1:3) - cos_jik*cdxp_ji(1:3))/r_ji
            dcos_jik_jk(1:3) = (cdxp_ji(1:3) - cos_jik*cdxp_jk(1:3))/r_jk
          end if

          do p = 1, g3_dim
            ang = (1.d0 + g3_lambda(p)*cos_jik)**g3_zeta(p)
            rad = dexp(-g3_eta(p)*(r2_ji + r2_jk + r2_ik))

            config_desc(iconf)%energy(p, ja) = config_desc(iconf)%energy(p, ja) + 0.5d0*2.d0**(1.d0 - g3_zeta(p))*ang*rad*fcut_ji*fcut_jk*fcut_ik

            if (desc_forces_local) then

              dang = g3_zeta(p)*(1.d0 + g3_lambda(p)*cos_jik)**(g3_zeta(p) - 1.d0)*g3_lambda(p)
              drad_ji(1:3) = -2.d0*g3_eta(p)*rad*dxp_ji(1:3)
              drad_jk(1:3) = -2.d0*g3_eta(p)*rad*dxp_jk(1:3)
              drad_ik(1:3) = -2.d0*g3_eta(p)*rad*dxp_ik(1:3)

              config_desc(iconf)%force(p, ja, ia_n, 1:3) = config_desc(iconf)%force(p, ja, ia_n, 1:3) + 0.5d0*2d0**(1d0 - g3_zeta(p))*( &
                                                           (dang*rad*dcos_jik_ji(1:3) + ang*(drad_ji(1:3) - drad_ik(1:3)))*fcut_ji*fcut_jk*fcut_ik + &
                                                           ang*rad*fcut_jk*(dfcut_ji(1:3)*fcut_ik - fcut_ji*dfcut_ik(1:3)))

              config_desc(iconf)%force(p, ja, ka_n, 1:3) = config_desc(iconf)%force(p, ja, ka_n, 1:3) + 0.5d0*2d0**(1d0 - g3_zeta(p))*( &
                                                           (dang*rad*dcos_jik_jk(1:3) + ang*(drad_jk(1:3) + drad_ik(1:3)))*fcut_ji*fcut_jk*fcut_ik + &
                                                           ang*rad*fcut_ji*(dfcut_jk(1:3)*fcut_ik + fcut_jk*dfcut_ik(1:3)))

            end if
            !config_desc(iconf)%force(p, ja, 0, 1:3) = config_desc(iconf)%force(p, ja, 0,  1:3) - 0.5d0*2d0**(1d0-g3_zeta(p))*( &
            !(dang*rad*(dcos_jik_ji(1:3) + dcos_jik_jk) + ang*(drad_ji(1:3) + drad_jk(1:3)) )*fcut_ji*fcut_jk*fcut_ik + &
            !   ang*rad*fcut_ik*(dfcut_ji(1:3)*fcut_jk + fcut_ji*dfcut_jk(1:3) )  )

          end do
        end do                  ! ka
        local_d_kind_neigh(ja, ia_n) = ia
      end do                  ! ia

      local_d_n_neigh(ja) = ia_n
    end do                  ! ja
    if (desc_forces_local) config_desc(iconf)%force(:, :, 0, :) = -SUM(config_desc(iconf)%force(:, :, 1:imm_neigh, :), dim=3)

    deallocate (xpnp)
  end subroutine compute_g3


  subroutine init_g3()
  use module_neigh_local, only: r_cut 
    use ml_in_ndm_module, only: rangml, g3_dim, strict_behler, n_g3_eta, n_g3_zeta, n_g3_lambda, &
                                g3_eta, g3_zeta, g3_lambda

    real(kind(0.d0)), dimension(:), allocatable  :: vd1
    real(kind(0.d0)), parameter      :: b2a = 0.52917721077d0
    integer  :: icount, p_e, p_l, p_z, ii

    g3_dim = n_g3_eta*n_g3_zeta*n_g3_lambda
    if (strict_behler) then
      g3_dim = 43
    end if

    if (allocated(g3_eta)) deallocate (g3_eta); allocate (g3_eta(g3_dim))
    if (allocated(g3_zeta)) deallocate (g3_zeta); allocate (g3_zeta(g3_dim))
    if (allocated(g3_lambda)) deallocate (g3_lambda); allocate (g3_lambda(g3_dim))


    if (strict_behler) then
      r_cut = 11.338d0*b2a

      if (allocated(vd1)) deallocate (vd1); allocate (vd1(g3_dim))
      vd1(:) = (/0.0001d0, 0.0001d0, 0.0001d0, 0.0001d0, &
                 0.0030d0, 0.0030d0, 0.0030d0, 0.0030d0, &
                 0.0080d0, 0.0080d0, 0.0080d0, 0.0080d0, &
                 0.0150d0, 0.0150d0, 0.0150d0, 0.0150d0, 0.0150d0, 0.0150d0, 0.0150d0, 0.0150d0, &
                 0.0250d0, 0.0250d0, 0.0250d0, 0.0250d0, 0.0250d0, 0.0250d0, 0.0250d0, 0.0250d0, &
                 0.0450d0, 0.0450d0, 0.0450d0, 0.0450d0, 0.0450d0, 0.0450d0, 0.0450d0, 0.0450d0, &
                 0.0800d0, 0.0800d0, 0.0800d0, 0.0800d0, 0.0800d0, 0.0800d0, 0.0800d0/)
      g3_eta(1:g3_dim) = vd1(1:g3_dim)/b2a**2

      do ii = 1, size(g3_lambda)
        g3_zeta(ii) = (-1.d0)**ii
      end do
      !g3_lambda(43)=1.d0
      g3_zeta(:) = (/1.000d0, 1.000d0, 2.000d0, 2.000d0, &
                     1.000d0, 1.000d0, 2.000d0, 2.000d0, &
                     1.000d0, 1.000d0, 2.000d0, 2.000d0, &
                     1.000d0, 1.000d0, 2.000d0, 2.000d0, 4.000d0, 4.000d0, 16.000d0, 16.000d0, &
                     1.000d0, 1.000d0, 2.000d0, 2.000d0, 4.000d0, 4.000d0, 16.000d0, 16.000d0, &
                     1.000d0, 1.000d0, 2.000d0, 2.000d0, 4.000d0, 4.000d0, 16.000d0, 16.000d0, &
                     1.000d0, 1.000d0, 2.000d0, 2.000d0, 4.000d0, 4.000d0, 16.000d0/)
    else
      icount = 0
      do p_e = 1, n_g3_eta
        do p_l = 1, n_g3_lambda
          do p_z = 1, n_g3_zeta
            icount = icount + 1
            g3_eta(icount) = 1.d-2 + dble(p_e - 1)*(0.80d0 - 1.d-2)/dble(n_g3_eta - 1)
            g3_lambda(icount) = -1.d0 + dble(p_l - 1)*2.d0
            g3_zeta(icount) = 2.d0**(p_z - 1)
          end do
        end do
      end do


      if (icount /= g3_dim) then
        if (rangml == 0) write (6, *) 'Error in init_g3'
        stop 'dimension of g3 is incorrect'
      end if
    end if
  end subroutine init_g3

end module
