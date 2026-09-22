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

module compute_afs_mod

  ! AFS computing the triangles centered on the atom j and the angle ijk.
  !
  !                 j--------i
  !                  \
  !                   \
  !                    \
  !                     k
  ! used quantities ji, kj and jik

#if(PARA)
  !use mpi
  use mld_mpi 
#endif

  USE module_kind_variables, ONLY: double

contains

  ! subroutine compute_afs(i_start_at, i_final_at, d_n_neigh, d_kind_neigh, local_afs_out, config_desc(iconf)%force, iconf)
  subroutine compute_afs(i_start_at, i_final_at, d_n_neigh, d_kind_neigh, iconf)

#ifdef MLD_NDM
    use gen_com_m, ONLY: A2cm, lperiod
    use gen_com_m_ml, ONLY: imm, bg, at
    use tab_imm_m_ml, ONLY: xp
#else
    use ondm_gen_com_m, ONLY: imm, A2cm, lperiod, bg, at, indi2
    use ondm_tab_imm_m, ONLY: iwmax2, xp
#endif
    use module_afs, only: n_rbf_afs, n_cheb, afs_dim, W_afs, coeff_rbf_afs, &
                          afs_type, afs_type_bartok
    use ml_in_ndm_module, ONLY: weighted, weighted_3ch, &
                                imm_neigh, desc_forces, linvisible
    use derived_types, only: config_real, config_desc
    use time_check_general, only: time, tot_time, debug_time, MY_MPI_WTIME
    use module_neigh_local, only: r_cut 
#ifdef MLD_NDM
    use notperiod_mod
#else
    use ondm_transform_coord, only: ondm_notperiod
#endif
    implicit none 

    integer, intent(in)  :: i_start_at, i_final_at
    integer, dimension(imm), intent(out)   :: d_n_neigh
    integer, dimension(imm, imm_neigh), intent(out)    :: d_kind_neigh
    !double precision,dimension(afs_dim, imm),intent(out)   :: config_desc(iconf)%energy
    !double precision,dimension(afs_dim, imm,0:imm_neigh, 3),intent(out) :: config_desc(iconf)%force
    integer, optional    :: iconf

    logical  :: small
    real(double), dimension(:, :), allocatable   :: xpnp
    real(double), dimension(3) :: dxp_ji, dxp_jk, dxp_ik, ds
    real(double), dimension(3) :: cdxp_ji, cdxp_jk, cdxp_ik
    real(double), dimension(n_rbf_afs)     :: phi_ji, dphi_ji, phi_jk, dphi_jk, rbf_ji, rbf_jk, drbf_ji, drbf_jk
    integer  :: ia, ja, ka, ia_n, ka_n, iw, i_desc, i_desc2, i_desc3
    integer  :: iw1, iw2
    integer  :: iz, p1, p2, p0, pmin, pmax
    double precision     :: r2_ji, r2_jk, r_ji, r_jk, r2_ik, r_ik, cos_kji, tmpr, term_local
    double precision, dimension(0:n_cheb)  :: chebT_kji, chebU_kji
    double precision, dimension(3)   :: dcos_kji_ji, dcos_kji_jk, dchebT_kji_ji, dchebT_kji_jk, term_local_3d
    double precision     :: factor_ia, factor_ia_3ch, factor_ja, factor_ja_3ch, factor_ka, factor_ka_3ch
    logical  :: desc_forces_local


    if ((i_start_at == 0) .and. (i_final_at == 0)) then
      d_n_neigh(:) = 0
      d_kind_neigh(:, :) = 0
      config_desc(iconf)%energy(:, :) = 0.d0
      !  config_desc(iconf)%force(:,:,:,:)=0.d0
      return
    end if
    desc_forces_local = desc_forces .and. (config_real(iconf)%has_force .or. config_real(iconf)%has_stress)
#ifndef MLD_NDM
    small = .false.
#endif
#ifdef MLD_NDM
!    if (present(iconf)) then
#else
    if (present(iconf)) then
#endif

      small = config_real(iconf)%small
#ifdef MLD_NDM
!    end if
!JPC xpnp ne sert a rien 
#else
    end if
#endif
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
!JPC
!    if (i_start_at == 1) iw2 = 0
!    if (i_start_at > 1) iw2 = iwmax2(i_start_at - 1)
    iw2=0
#else
    if (i_start_at == 1) iw2 = 0
    if (i_start_at > 1) iw2 = iwmax2(i_start_at - 1)
#endif
    if (debug_time) time(1) = MY_MPI_WTIME()

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
      !end   small box or not 1/

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
        phi_ji(:) = 0.d0
        ia_n = ia_n + 1

        do p1 = 1, n_rbf_afs
          !phi_ji(p1) = dsqrt((2.d0*p1+5.d0)*r_cut**(-2.d0*p1-5.d0))*(r_cut-r_ji)**(p1+2.d0)
          tmpr = (r_cut - r_ji)**(p1 + 1)
          phi_ji(p1) = coeff_rbf_afs(p1)*tmpr*(r_cut - r_ji)
          !phi_ji(p1) = coeff_rbf_afs(p1)*(r_cut-r_ji)**(p1+2.d0)
          !if (desc_forces_local) dphi_ji(p1)=-coeff_rbf_afs(p1)*(p1+2.d0)*(r_cut-r_ji)**(p1+1.d0)
          if (desc_forces_local) dphi_ji(p1) = -coeff_rbf_afs(p1)*(p1 + 2)*tmpr
        end do
        rbf_ji(:) = matmul(W_afs(:, :), phi_ji(:))
        if (desc_forces_local) drbf_ji(:) = matmul(W_afs(:, :), dphi_ji(:))
        cdxp_ji(:) = dxp_ji(:)/r_ji


        if (debug_time) time(3) = MY_MPI_WTIME()

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
            if (weighted_3ch) factor_ka_3ch = config_real(iconf)%weight_per_type_3ch(config_real(iconf)%itype(ka))
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
          !why_not if (r_ik >= r_cut) cycle
          if (r_ik <= 1d-20) cycle

          do p1 = 1, n_rbf_afs
            tmpr = (r_cut - r_jk)**(p1 + 1)
            !phi_jk(p1) = coeff_rbf_afs(p1)*(r_cut-r_jk)**(p1+2.d0)
            phi_jk(p1) = coeff_rbf_afs(p1)*(r_cut - r_jk)*tmpr
            !if (desc_forces) dphi_jk(p1)=-dsqrt((2.d0*p1+5.d0)*r_cut**(-2.d0*p1-5.d0))*(p1+2.d0)*(r_cut-r_jk)**(p1+1.d0)
            if (desc_forces) dphi_jk(p1) = -coeff_rbf_afs(p1)*dble(p1 + 2)*tmpr
          end do
          rbf_jk(:) = matmul(W_afs(:, :), phi_jk(:))
          if (desc_forces_local) drbf_jk(:) = matmul(W_afs(:, :), dphi_jk(:))

          cdxp_jk(:) = dxp_jk(:)/r_jk
          cdxp_ik(:) = dxp_ik(:)/r_ik
          cos_kji = dot_product(cdxp_ji(1:3), cdxp_jk(1:3))

          if (desc_forces_local) dcos_kji_ji(1:3) = (cdxp_jk(1:3) - cos_kji*cdxp_ji(1:3))/r_ji
          if (desc_forces_local) dcos_kji_jk(1:3) = (cdxp_ji(1:3) - cos_kji*cdxp_jk(1:3))/r_jk

          chebT_kji(0) = 1.d0
          chebT_kji(1) = cos_kji
          chebU_kji(0) = 1.d0
          chebU_kji(1) = 2.d0*cos_kji
          if (n_cheb .ge. 2) then
            do p2 = 2, n_cheb
              chebT_kji(p2) = 2.d0*cos_kji*chebT_kji(p2 - 1) - chebT_kji(p2 - 2)
              chebU_kji(p2) = 2.d0*cos_kji*chebU_kji(p2 - 1) - chebU_kji(p2 - 2)
            end do
          end if

          i_desc = 0
          do p0 = 1, n_rbf_afs
            if (afs_type == afs_type_bartok) then
              pmin = p0
              pmax = p0
            else
              pmin = 1
              pmax = n_rbf_afs
            end if

            do p1 = pmin, pmax
              do p2 = 0, n_cheb
                if (desc_forces_local) then 
                  if (p2 == 0) then
                    dchebT_kji_ji(1:3) = 0.d0
                    dchebT_kji_jk(1:3) = 0.d0
                  else
                    dchebT_kji_ji(1:3) = p2*chebU_kji(p2 - 1)*dcos_kji_ji(1:3)
                    dchebT_kji_jk(1:3) = p2*chebU_kji(p2 - 1)*dcos_kji_jk(1:3)
                  end if
                end if 

                i_desc = i_desc + 1
                !oldWesley config_desc(iconf)%energy(i_desc,ja) = config_desc(iconf)%energy(i_desc,ja) + factor_ia*factor_ka* rbf_ji(p1,ia,ja) * rbf(p1,ka,ja) * chebT_kji(p2)
                tmpr = rbf_ji(p0)*rbf_jk(p1)
                !term_local=rbf_ji(p0) * rbf_jk(p1) * chebT_kji(p2)/2.d0
                term_local = tmpr*chebT_kji(p2)/2.d0
                config_desc(iconf)%energy(i_desc, ja) = config_desc(iconf)%energy(i_desc, ja) + term_local
                if (weighted) then
                  i_desc2 = i_desc + afs_dim
                  i_desc3 = i_desc2 + afs_dim
                  config_desc(iconf)%energy(i_desc2, ja) = config_desc(iconf)%energy(i_desc2, ja) + factor_ia*factor_ka*term_local
                  if (weighted_3ch) config_desc(iconf)%energy(i_desc3, ja) = config_desc(iconf)%energy(i_desc3, ja) + factor_ia_3ch*factor_ka_3ch*term_local
                end if

                if (desc_forces_local) then
                  term_local_3d(1:3) = (rbf_ji(p0)*rbf_jk(p1)*dchebT_kji_ji(1:3) + drbf_ji(p0)*cdxp_ji(1:3)*rbf_jk(p1)*chebT_kji(p2))         ! /2.d0
                  !term_local_3d(1:3)= rbf_ji(p0)*rbf_jk(p1)*(dchebT_kji_ji(1:3) + dchebT_kji_jk(1:3))  + &
                  !                    (drbf_ji(p0)*cdxp_ji(1:3)*rbf_jk(p1) + rbf_ji(p0)*cdxp_jk(1:3)*drbf_jk(p1)) *chebT_kji(p2)
                  !term_local_3d(1:3)= tmpr*(dchebT_kji_ji(1:3) + dchebT_kji_jk(1:3))  + &
                  !                    (drbf_ji(p0)*cdxp_ji(1:3)*rbf_jk(p1) + rbf_ji(p0)*cdxp_jk(1:3)*drbf_jk(p1)) *chebT_kji(p2)
                  config_desc(iconf)%force(i_desc, ja, ia_n, 1:3) = config_desc(iconf)%force(i_desc, ja, ia_n, 1:3) + term_local_3d(1:3)
                  if (weighted) then
                    config_desc(iconf)%force(i_desc2, ja, ia_n, 1:3) = config_desc(iconf)%force(i_desc2, ja, ia_n, 1:3) + factor_ia*factor_ka*term_local_3d(1:3)
                    if (weighted_3ch) config_desc(iconf)%force(i_desc3, ja, ia_n, 1:3) = config_desc(iconf)%force(i_desc3, ja, ia_n, 1:3) + factor_ia_3ch*factor_ka_3ch*term_local_3d(1:3)
                  end if
                end if

              end do                  ! p2 - cos
            end do                  ! p1
          end do                  ! p0
        end do                  ! ka of the  neigh of ja

        if (debug_time) then
          time(4) = MY_MPI_WTIME()
          tot_time(2) = tot_time(2) + time(4) - time(3)
        end if

        d_kind_neigh(ja, ia_n) = ia
        if (desc_forces_local) config_desc(iconf)%force(:, ja, 0, :) = config_desc(iconf)%force(:, ja, 0, :) - config_desc(iconf)%force(:, ja, ia_n, :)      ! /0.5d0
      end do                  ! ia_n loop over neigh of ja

      if (weighted) then
        config_desc(iconf)%energy(1 + afs_dim:2*afs_dim, ja) = config_desc(iconf)%energy(1 + afs_dim:2*afs_dim, ja)*factor_ja
        if (weighted_3ch) config_desc(iconf)%energy(2*afs_dim + 1:3*afs_dim, ja) = config_desc(iconf)%energy(2*afs_dim + 1:3*afs_dim, ja)*factor_ja_3ch
      end if
      if (desc_forces_local) then
        if (weighted) then
          config_desc(iconf)%force(1 + afs_dim:2*afs_dim, ja, :, :) = config_desc(iconf)%force(1 + afs_dim:2*afs_dim, ja, :, :)*factor_ja
          if (weighted_3ch) config_desc(iconf)%force(2*afs_dim + 1:3*afs_dim, ja, :, :) = config_desc(iconf)%force(2*afs_dim + 1:3*afs_dim, ja, :, :)*factor_ja_3ch
        end if
      end if
      d_n_neigh(ja) = ia_n
    end do                  ! ja main loop
    deallocate (xpnp)

    if (debug_time) then
      time(2) = MY_MPI_WTIME()
      tot_time(1) = tot_time(1) + time(2) - time(1)
    end if

    if (debug_time) then
      call repport_time(2, 0.d0, tot_time(1), "ML: full")
      call repport_time(4, 0.d0, tot_time(2), "ML: loop on k ")
    end if

  end subroutine compute_afs



  subroutine init_afs_rbf()

    use module_neigh_local, only: r_cut 
    use module_afs, only: n_rbf_afs, n_cheb, W_afs, afs_type, afs_type_bartok, afs_type_homemade, &
                          afs_dim, coeff_rbf_afs


    real(kind(0.d0)), dimension(n_rbf_afs, n_rbf_afs)  :: S, V
    real(kind(0.d0)), dimension(n_rbf_afs) :: L
    !local
    integer  :: p, q, nb, ilaenv, lwork


    if (allocated(W_afs)) deallocate (W_afs)
    allocate (W_afs(n_rbf_afs, n_rbf_afs))
    if (allocated(coeff_rbf_afs)) deallocate (coeff_rbf_afs)
    allocate (coeff_rbf_afs(n_rbf_afs))


    if (afs_type == afs_type_bartok) afs_dim = n_rbf_afs*(n_cheb + 1)
    if (afs_type == afs_type_homemade) afs_dim = n_rbf_afs**2*(n_cheb + 1)

    do q = 1, n_rbf_afs
      do p = 1, n_rbf_afs
        S(p, q) = dsqrt((2.d0*dble(p) + 5.d0)*(2.d0*dble(q) + 5.d0))/(dble(p + q) + 5.d0)
      end do
    end do

    nb = ilaenv(1, 'DSYTRD', 'L', n_rbf_afs, -1, -1, -1)
    lwork = (nb + 2)*n_rbf_afs
    call diagsym(S, n_rbf_afs, lwork, L)
    V(:, :) = 0d0
    do p = 1, n_rbf_afs
      if (L(p) == 0.d0) then
        !if (abs(L(p) - 0.d0).lt.1.d-15)  then
        V(p, p) = 0.d0
      else
        V(p, p) = L(p)/dabs(L(p))*dabs(L(p))**(-0.5)
      end if

      coeff_rbf_afs(p) = dsqrt((2.d0*p + 5.d0)*r_cut**(-2.d0*p - 5.d0))
    end do

    W_afs(:, :) = matmul(S(:, :), matmul(V(:, :), transpose(S(:, :))))

  end subroutine init_afs_rbf


  subroutine diagsym(A, n_rbf, lwork, L)
    ! input - the symmetric matrix A
    ! output - the ortho-normalized vectors that diaonalize A and the eigenvalues L
    use ml_in_ndm_module, only: rangml

    integer, intent(in)  :: n_rbf, lwork
    double precision, dimension(n_rbf, n_rbf), intent(inout) :: A
    double precision, dimension(n_rbf), intent(out)    :: L
    double precision, dimension(lwork)     :: work
    integer  :: info

    call dsyev('V', 'L', n_rbf, A, n_rbf, L, work, lwork, info)

    if (rangml == 0) then
      if (.not. (info == 0)) then
        write (6, *) 'ML: WARNING the nomalization of the overlap matrix is WRONG in diagsym compute_afs.F90'
      end if
    end if

  end subroutine diagsym

end module
