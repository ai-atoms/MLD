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

module compute_pow_so3_mod
#if(PARA)
  !use mpi
  use mld_mpi 
#endif
  USE module_kind_variables, ONLY: double, kind_double, kind_double_complex
  implicit none
contains

  ! subroutine compute_pow_so3(i_start_at,i_final_at,d_n_neigh, d_kind_neigh, local_pow_so3_out,local_pow_so3_deriv_out, iconf)
  subroutine compute_pow_so3(i_start_at, i_final_at, d_n_neigh, d_kind_neigh, iconf)

#ifdef MLD_NDM
    use gen_com_m, ONLY: A2cm, lperiod
    use gen_com_m_ml, ONLY: imm, bg, at
    use tab_imm_m_ml, ONLY: xp
#else
    use ondm_gen_com_m, ONLY: imm, A2cm, lperiod, bg, at, indi2
    use ondm_tab_imm_m, ONLY: iwmax2, xp
#endif
    use angular_functions, only: spherical_harm, grad_spherical_harm
    use derived_types, only: config_real, config_desc
    use ml_in_ndm_module, ONLY: imm_neigh, one_pi,  l_max, &
                                desc_forces, linvisible, weighted, weighted_3ch
    use module_neigh_local, only: r_cut 
    use time_check_general, only: time, tot_time, debug_time, MY_MPI_WTIME
    use module_so3, only: n_rbf_so3, pow_so3_dim, clmn, dclmn, ini_rbf_so3, &
                          clmn_w, dclmn_w, clmn_w_3ch, dclmn_w_3ch

#ifdef MLD_NDM
  use notperiod_mod
#else
    use ondm_transform_coord, only: ondm_notperiod
#endif
    integer, intent(in)  :: i_start_at, i_final_at
    integer, dimension(imm), intent(out)   :: d_n_neigh
    integer, dimension(imm, imm_neigh), intent(out)    :: d_kind_neigh
    integer, optional    :: iconf

    logical  :: small
    real(double), dimension(:, :), allocatable   :: xpnp
    real(double), dimension(3) :: dxp_ji, ds
    integer  :: ia, ja, iw, iw1, iw2, ia_n, i_count

    integer  :: p, l, m
    real(double)   :: r2_ji, r_ji, tmp_fact

    real(double), dimension(pow_so3_dim)   :: tmp_pow_so3_out
    real(double), dimension(pow_so3_dim, 0:imm_neigh, 3)     :: tmp_pow_so3_deriv_out

    real(double), dimension(pow_so3_dim)   :: tmp_pow_so3_out_w
    real(double), dimension(pow_so3_dim, 0:imm_neigh, 3)     :: tmp_pow_so3_deriv_out_w

    real(double), dimension(pow_so3_dim)   :: tmp_pow_so3_out_w_3ch
    real(double), dimension(pow_so3_dim, 0:imm_neigh, 3)     :: tmp_pow_so3_deriv_out_w_3ch

    real(double), dimension(n_rbf_so3)     :: phi_ji
    real(double)   :: factor_ia, factor_ia_3ch, factor_ja, factor_ja_3ch
    logical  :: desc_forces_local

    if ((i_start_at == 0) .and. (i_final_at == 0)) then
      d_n_neigh(:) = 0
      d_kind_neigh(:, :) = 0
      config_desc(iconf)%energy(:, :) = 0.d0
      ! config_desc(iconf)%force(:,:,:,:) = 0.d0
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
        clmn_w(:, :) = (0.d0, 0.d0)
        tmp_pow_so3_out_w(:) = 0.d0
        tmp_pow_so3_deriv_out_w(:, :, :) = 0.d0
        if (weighted_3ch) then
          clmn_w_3ch(:, :) = (0.d0, 0.d0)
          tmp_pow_so3_out_w_3ch(:) = 0.d0
          tmp_pow_so3_deriv_out_w_3ch(:, :, :) = 0.d0
        end if
      end if

      ia_n = 0

      if (debug_time) then
        time(4) = MY_MPI_WTIME()
        tot_time(2) = tot_time(2) + time(4) - time(3)
      end if


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
      ! end   small box or not 1/

      do iw = iw1, iw2
        !begin small box or not 2/
#if(PARA)
        if (debug_time) time(7) = MY_MPI_WTIME()
#else
        if (debug_time) time(7) = 0                      ! TODO
#endif
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

        factor_ia = 1.d0
        factor_ia_3ch = 1.d0
        if (weighted) then
          if (linvisible) then
            if (config_real(iconf)%invisible_per_type(config_real(iconf)%itype(ia))) cycle
          end if
          factor_ia = config_real(iconf)%weight_per_type(config_real(iconf)%itype(ia))
          if (weighted_3ch) factor_ia_3ch = config_real(iconf)%weight_per_type_3ch(config_real(iconf)%itype(ia))
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
!!$          dxp_ji(:) = MatMul(at(:, :), ds(:))/A2cm
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
          dxp_ji(:) = MatMul(at(:, :), ds(:))/A2cm
          r2_ji = Sum(dxp_ji(1:3)**2)
          r_ji = dsqrt(r2_ji)
        end if
#endif

        phi_ji(:) = 0.d0
        if (r_ji >= r_cut) cycle
        ia_n = ia_n + 1

        if (debug_time) then
          time(8) = MY_MPI_WTIME()
          tot_time(5) = tot_time(5) + time(8) - time(7)
        end if

        ! this gives clmn and d_clmn ...
        call spherical_3d(ia_n, factor_ia, factor_ia_3ch, desc_forces_local, dxp_ji, r_ji)

        d_kind_neigh(ja, ia_n) = ia

        if (debug_time) then
          time(9) = MY_MPI_WTIME()
          tot_time(6) = tot_time(6) + time(9) - time(8)
        end if

      end do                  ! ia end of neighbours iterations ...

      ! update the componenet from the central atom ....
      i_count = 0
      do p = ini_rbf_so3, n_rbf_so3
        do l = 0, l_max
          i_count = i_count + 1
          tmp_fact =  dsqrt( (2.d0*dble(l)+1) / (4.d0 * one_pi) )
          clmn(0, i_count) = clmn(0, i_count) + tmp_fact
          if (weighted) then
            clmn_w(0, i_count) = clmn_w(0, i_count) + tmp_fact*factor_ja
            if (weighted_3ch) clmn_w_3ch(0, i_count) = clmn_w_3ch(0, i_count) + tmp_fact*factor_ja_3ch
          end if
        end do
      end do

      d_n_neigh(ja) = ia_n
      if (debug_time) then
        time(5) = MY_MPI_WTIME()
        tot_time(3) = tot_time(3) + time(5) - time(4)
      end if
      i_count = 0
      do p = ini_rbf_so3, n_rbf_so3
        do l = 0, l_max
          i_count = i_count + 1
          do m = -l, l

            tmp_pow_so3_out(i_count) = tmp_pow_so3_out(i_count) + real(dconjg(clmn(m, i_count))*clmn(m, i_count), kind(0.d0))
            if (weighted) then
              tmp_pow_so3_out_w(i_count) = tmp_pow_so3_out_w(i_count) + real(dconjg(clmn_w(m, i_count))*clmn_w(m, i_count), kind(0.d0))
              if (weighted_3ch) then
                tmp_pow_so3_out_w_3ch(i_count) = tmp_pow_so3_out_w_3ch(i_count) + real(dconjg(clmn_w_3ch(m, i_count))*clmn_w_3ch(m, i_count), kind(0.d0))
              end if
            end if

            do ia = 1, d_n_neigh(ja)
              if (desc_forces_local) tmp_pow_so3_deriv_out(i_count, ia, 1:3) = tmp_pow_so3_deriv_out(i_count, ia, 1:3) + &
                                                                               2.d0*real(dconjg(clmn(m, i_count))*dclmn(m, i_count, ia, 1:3), kind(0.d0))
            end do
            if (weighted) then
              do ia = 1, d_n_neigh(ja)
                if (desc_forces_local) tmp_pow_so3_deriv_out_w(i_count, ia, 1:3) = tmp_pow_so3_deriv_out_w(i_count, ia, 1:3) + &
                                                                               2.d0*real(dconjg(clmn_w(m, i_count))*dclmn_w(m, i_count, ia, 1:3), kind(0.d0))
              end do
              if (weighted_3ch) then
                do ia = 1, d_n_neigh(ja)
                  if (desc_forces_local) tmp_pow_so3_deriv_out_w_3ch(i_count, ia, 1:3) = tmp_pow_so3_deriv_out_w_3ch(i_count, ia, 1:3) + &
                                                                               2.d0*real(dconjg(clmn_w_3ch(m, i_count))*dclmn_w_3ch(m, i_count, ia, 1:3), kind(0.d0))
                end do
              end if
            end if
          end do                  ! m
        end do                  ! p
      end do                  ! l
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
  end subroutine compute_pow_so3


  subroutine init_pow_so3()

    use ml_in_ndm_module, only: l_max, imm_neigh, &
                                weighted, weighted_3ch
    use module_so3, only: n_rbf_so3, pow_so3_dim,  radial_pow_so3, radial_bartok, radial_sgg, clmn, dclmn, &
                          clmn_w, dclmn_w, clmn_w_3ch, dclmn_w_3ch,  ini_rbf_so3


    if (radial_pow_so3 == radial_bartok) then
      ini_rbf_so3 = 1
      pow_so3_dim = int((1 + l_max))*n_rbf_so3
    end if
    if (radial_pow_so3 == radial_sgg) then
      ini_rbf_so3 = 0
      pow_so3_dim = int((1 + l_max))*(n_rbf_so3 + 1)
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

  end subroutine init_pow_so3


  subroutine init_pow_so3_rbf()
    use module_neigh_local, only: r_cut 
    use module_so3, only: n_rbf_so3, coeff_rbf_so3, W_pow_so3, &
                          chebT, chebU, d_chebT, &
                          scgg, d_scgg, radial_pow_so3, radial_bartok, &
                          radial_pow_so3, radial_sgg
    use compute_afs_mod


    real(kind(0.d0)), dimension(n_rbf_so3, n_rbf_so3)  :: S, V
    real(kind(0.d0)), dimension(n_rbf_so3) :: L
    !local
    integer  :: p, q, nb, ilaenv, lwork


    if (radial_pow_so3 == radial_bartok) then
      if (allocated(W_pow_so3)) deallocate (W_pow_so3); allocate (W_pow_so3(n_rbf_so3, n_rbf_so3))
      if (allocated(coeff_rbf_so3)) deallocate (coeff_rbf_so3); allocate (coeff_rbf_so3(n_rbf_so3))

      do q = 1, n_rbf_so3
        do p = 1, n_rbf_so3
          !Bartok paper
          S(p, q) = dsqrt((2.d0*dble(p) + 5.d0)*(2.d0*dble(q) + 5.d0))/(dble(p + q) + 5.d0)
          !trail version:
          ! S(p,q)=dsqrt((2.d0*dble(p)+5.d0)*(2.d0*dble(p)+6.d0)*(2.d0*dble(p)+7.d0)*  &
          !              (2.d0*dble(q)+5.d0)*(2.d0*dble(q)+6.d0)*(2.d0*dble(q)+7.d0) ) &
          !              /                                                             &
          !              ( (5.d0+dble(p)+dble(q)) * (6.d0+dble(p)+dble(q)) * (7.d0+dble(p)+dble(q)) )
        end do
      end do

      nb = ilaenv(1, 'DSYTRD', 'L', n_rbf_so3, -1, -1, -1)
      lwork = (nb + 2)*n_rbf_so3
      call diagsym(S, n_rbf_so3, lwork, L)
      V(:, :) = 0d0
      do p = 1, n_rbf_so3
        if (L(p) == 0.d0) then
          V(p, p) = 0.d0
        else
          !V(p, p) = L(p)/dabs(L(p))*dabs(L(p))**(-0.5)
          V(p, p) = sign(1.d0,L(p))*dabs(L(p))**(-0.5)
        end if
        ! Bartok paper
        coeff_rbf_so3(p) = dsqrt((2.d0*p + 5.d0)*r_cut**(-2.d0*p - 5.d0))
        !trial version:
        !coeff_rbf_so3(p) =  dsqrt( 0.5d0 * (2.d0*p+5.d0) * (2.d0*p+6.d0) * (2.d0*p+7.d0) * r_cut**(-2*p-7))

      end do

      W_pow_so3(:, :) = matmul(S(:, :), matmul(V(:, :), transpose(S(:, :))))
    end if

    if (radial_pow_so3 == radial_sgg) then
      if (allocated(chebT)) deallocate (chebT); allocate (chebT(0:n_rbf_so3))
      if (allocated(d_chebT)) deallocate (d_chebT); allocate (d_chebT(0:n_rbf_so3))
      if (allocated(chebU)) deallocate (chebU); allocate (chebU(0:n_rbf_so3))

      if (allocated(scgg)) deallocate (scgg); allocate (scgg(0:n_rbf_so3))
      if (allocated(d_scgg)) deallocate (d_scgg); allocate (d_scgg(0:n_rbf_so3))
    end if

  end subroutine init_pow_so3_rbf



  subroutine spherical_3d(ia_n, factor_ia, factor_ia_3ch, desc_forces_local, dxp, rr)

    use ml_in_ndm_module, only: l_max,   &
                                weighted, weighted_3ch
    use module_so3, only: n_rbf_so3, dclmn, clmn, coeff_rbf_so3, W_pow_so3, radial_pow_so3, &
                          radial_bartok, radial_sgg, radial_ace_bessel, &
                          scgg, d_scgg, &
                          ini_rbf_so3, end_rbf_so3, &
                          dclmn_w, clmn_w, dclmn_w_3ch, clmn_w_3ch, &
                          radial_ace, d_radial_ace
    use module_neigh_local, only: r_cut 
    use mod_radial_functions, only:  paftouny_so3, sc_stefanodg
    use angular_functions, only: spherical_harm, grad_spherical_harm
    use time_check_general, only: debug_time,  tot_time, debug_time, MY_MPI_WTIME
    implicit none
    real(kind_double), intent(in)    :: rr, factor_ia, factor_ia_3ch
    real(kind_double), dimension(3), intent(in)  :: dxp
    integer, intent(in)  :: ia_n
    logical, intent(in)  :: desc_forces_local
    !complex(kind_double_complex), dimension(-l_max:l_max,pow_so3_dim) :: cmm
    real(kind_double), dimension(:), allocatable :: phi_ji, dphi_ji, rbf_ji, drbf_ji
    real(kind_double)    :: tmp_rr, tmp_rrp, tmp_rbf, tmp_rbf_w, tmp_rbf_w_3ch, &
                            dtmp_rbf, dtmp_rbf_w, dtmp_rbf_w_3ch, tmp_cos(3)
    !------------
    real(kind_double), parameter     :: alpha_gemv = 1.d0, beta_gemv = 0.d0
    integer, parameter   :: incx_gemv = 1, incy_gemv = 1
    !-------------
    integer  :: p, l, i_count, m
    real(kind_double)    :: time0, time1, time2, time3
    double complex, dimension(0:l_max, -l_max:l_max)   :: sharm
    double complex, dimension(0:l_max, -l_max:l_max, 3)      :: grad_sharm

    if (debug_time) time0 = MY_MPI_WTIME()
    if (allocated(phi_ji)) deallocate (phi_ji)
    allocate (phi_ji(n_rbf_so3))
    if (allocated(dphi_ji)) deallocate (dphi_ji)
    allocate (dphi_ji(n_rbf_so3))
    if (allocated(rbf_ji)) deallocate (rbf_ji)
    allocate (rbf_ji(n_rbf_so3))
    if (allocated(drbf_ji)) deallocate (drbf_ji)
    allocate (drbf_ji(n_rbf_so3))

    tmp_cos(1:3) = dxp(1:3)/rr

    if (radial_pow_so3 == radial_bartok) then
      rbf_ji(:)  = 0.d0
      drbf_ji(:) = 0.d0
      tmp_rr = r_cut - rr
      do p = 1, n_rbf_so3
        tmp_rrp = tmp_rr**(p + 1)
        phi_ji(p) = coeff_rbf_so3(p)*tmp_rr*tmp_rrp
        if (desc_forces_local) dphi_ji(p) = -coeff_rbf_so3(p)*(p + 2)*tmp_rrp
      end do
      ! y = alpha*A*x + beta * y
      !call dgemv(trans, m, n, alpha, a, lda, x, incx, beta, y, incy)

      !rbf_ji(:) =matmul(W_pow_so3(:,:) , phi_ji(:))
      ! y = \alpha *W phi_ji + \beta y
      call dgemv('N', n_rbf_so3, n_rbf_so3, alpha_gemv, W_pow_so3, n_rbf_so3, phi_ji, incx_gemv, beta_gemv, rbf_ji, incy_gemv)
      if (desc_forces_local) then
        !drbf_ji(:)=matmul(W_pow_so3(:,:) ,dphi_ji(:))
        call dgemv('N', n_rbf_so3, n_rbf_so3, alpha_gemv, W_pow_so3, n_rbf_so3, dphi_ji, incx_gemv, beta_gemv, drbf_ji, incy_gemv)
      end if
    end if

    if (radial_pow_so3 == radial_sgg) then
      tmp_rr = 2.d0*rr/r_cut - 1.d0
      call paftouny_so3(desc_forces_local, n_rbf_so3, tmp_rr)
      call sc_stefanodg(desc_forces_local, n_rbf_so3, rr)
      if (allocated(rbf_ji)) deallocate (rbf_ji); allocate (rbf_ji(0:n_rbf_so3))
      if (allocated(drbf_ji)) deallocate (drbf_ji); allocate (drbf_ji(0:n_rbf_so3))
    end if


    if (debug_time) then
      time1 = MY_MPI_WTIME()
      tot_time(8) = tot_time(8) + time1 - time0
    end if

    ! compute spherical functions
    do l = 0, l_max
      do m = -l, l
        sharm(l, m) = spherical_harm(l, m, dxp(:))
        if (desc_forces_local) grad_sharm(l, m, :) = grad_spherical_harm(l, m, dxp(:))
      end do
    end do

    if (debug_time) then
      time2 = MY_MPI_WTIME()
      tot_time(9) = tot_time(9) + time2 - time1
    end if

    i_count = 0
    if (radial_pow_so3 == radial_bartok) then
    do p = ini_rbf_so3, n_rbf_so3
      tmp_rbf = rbf_ji(p)
      tmp_rbf_w = factor_ia * tmp_rbf
      tmp_rbf_w_3ch = factor_ia_3ch * tmp_rbf
      dtmp_rbf = drbf_ji(p)
      dtmp_rbf_w = factor_ia * dtmp_rbf
      dtmp_rbf_w_3ch = factor_ia_3ch * dtmp_rbf
      do l = 0, l_max
        i_count = i_count + 1
        do m = -l, l
          !clmn(m, i_count) = clmn(m, i_count) + rbf_ji(p)*sharm(l, m)
          clmn(m, i_count) = clmn(m, i_count) + tmp_rbf*sharm(l, m)
          !if (desc_forces_local) dclmn(m, i_count, ia_n, 1:3) = drbf_ji(p)*tmp_cos(1:3)*sharm(l, m) + &
          !                                                      rbf_ji(p)*grad_sharm(l, m, 1:3)
          if (desc_forces_local) dclmn(m, i_count, ia_n, 1:3) = dtmp_rbf*tmp_cos(1:3)*sharm(l, m) + &
                                                                tmp_rbf*grad_sharm(l, m, 1:3)
          if (weighted) then
            clmn_w(m, i_count) = clmn_w(m, i_count) + tmp_rbf_w*sharm(l, m)
            if (desc_forces_local) dclmn_w(m, i_count, ia_n, 1:3) = dtmp_rbf_w*tmp_cos(1:3)*sharm(l, m) + &
                                                                tmp_rbf_w*grad_sharm(l, m, 1:3)
            if (weighted_3ch) then
              clmn_w_3ch(m, i_count) = clmn_w_3ch(m, i_count) + tmp_rbf_w_3ch*sharm(l, m)
              if (desc_forces_local) dclmn_w_3ch(m, i_count, ia_n, 1:3) = dtmp_rbf_w_3ch*tmp_cos(1:3)*sharm(l, m) + &
                                                                tmp_rbf_w_3ch*grad_sharm(l, m, 1:3)

            end if
          end if
        end do                  ! m
      end do                  ! l
    end do                  ! p
    end if


    if (radial_pow_so3 == radial_sgg) then
    do p = ini_rbf_so3, n_rbf_so3
      tmp_rbf = scgg(p)
      tmp_rbf_w = factor_ia * tmp_rbf
      tmp_rbf_w_3ch = factor_ia_3ch * tmp_rbf
      dtmp_rbf = d_scgg(p)
      dtmp_rbf_w = factor_ia * dtmp_rbf
      dtmp_rbf_w_3ch = factor_ia_3ch * dtmp_rbf
      do l = 0, l_max
        i_count = i_count + 1
        do m = -l, l
          clmn(m, i_count) = clmn(m, i_count) + tmp_rbf*sharm(l, m)
          if (desc_forces_local) dclmn(m, i_count, ia_n, 1:3) = dtmp_rbf*tmp_cos(1:3)*sharm(l, m) + &
                                                                tmp_rbf*grad_sharm(l, m, 1:3)
        end do
        if (weighted) then
          do m = -l, l
            clmn_w(m, i_count) = clmn_w(m, i_count) + tmp_rbf_w*sharm(l, m)
            if (desc_forces_local) dclmn_w(m, i_count, ia_n, 1:3) = dtmp_rbf_w*tmp_cos(1:3)*sharm(l, m) + &
                                                                tmp_rbf_w*grad_sharm(l, m, 1:3)
          end do
          if (weighted_3ch) then
            do m = -l, l
              clmn_w_3ch(m, i_count) = clmn_w_3ch(m, i_count) + tmp_rbf_w_3ch*sharm(l, m)
              if (desc_forces_local) dclmn_w_3ch(m, i_count, ia_n, 1:3) = dtmp_rbf_w_3ch*tmp_cos(1:3)*sharm(l, m) + &
                                                                tmp_rbf_w_3ch*grad_sharm(l, m, 1:3)
            end do
          end if
        end if
      end do                  ! l
    end do                  ! p
    end if


    if (radial_pow_so3 == radial_ace_bessel) then
    do p = ini_rbf_so3, end_rbf_so3
      tmp_rbf = radial_ace(p)
      tmp_rbf_w = factor_ia * tmp_rbf
      tmp_rbf_w_3ch = factor_ia_3ch * tmp_rbf
      dtmp_rbf = d_radial_ace(p)
      dtmp_rbf_w = factor_ia * dtmp_rbf
      dtmp_rbf_w_3ch = factor_ia_3ch * dtmp_rbf
      do l = 0, l_max
        i_count = i_count + 1
        do m = -l, l
          clmn(m, i_count) = clmn(m, i_count) + tmp_rbf*sharm(l, m)
          if (desc_forces_local) dclmn(m, i_count, ia_n, 1:3) = dtmp_rbf*tmp_cos(1:3)*sharm(l, m) + &
                                                                tmp_rbf*grad_sharm(l, m, 1:3)
        end do
        if (weighted) then
          do m = -l, l
            clmn_w(m, i_count) = clmn_w(m, i_count) + tmp_rbf_w*sharm(l, m)
            if (desc_forces_local) dclmn_w(m, i_count, ia_n, 1:3) = dtmp_rbf_w*tmp_cos(1:3)*sharm(l, m) + &
                                                                tmp_rbf_w*grad_sharm(l, m, 1:3)
          end do
          if (weighted_3ch) then
            do m = -l, l
              clmn_w_3ch(m, i_count) = clmn_w_3ch(m, i_count) + tmp_rbf_w_3ch*sharm(l, m)
              if (desc_forces_local) dclmn_w_3ch(m, i_count, ia_n, 1:3) = dtmp_rbf_w_3ch*tmp_cos(1:3)*sharm(l, m) + &
                                                                tmp_rbf_w_3ch*grad_sharm(l, m, 1:3)
            end do
          end if
        end if
      end do                  ! l
    end do                  ! p
    end if




    if (debug_time) then
      time3 = MY_MPI_WTIME()
      tot_time(10) = tot_time(10) + time3 - time2
    end if
  end subroutine spherical_3d




end module
