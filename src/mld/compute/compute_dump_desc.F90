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

module compute_dump_desc_mod
  use module_kind_variables, ONLY: double

  implicit none

contains

  ! subroutine compute_bispectrum_so4(i_start_at,i_final_at,d_n_neigh, d_kind_neigh, local_bispectrum_so4_out,local_bispectrum_so4_deriv_out, iconf)
  subroutine compute_dump_desc(i_start_at, i_final_at, d_n_neigh, d_kind_neigh, iconf)

#ifdef MLD_NDM
    use gen_com_m, ONLY: A2cm, lperiod
    use gen_com_m_ml, ONLY: imm, bg, at
    use tab_imm_m_ml, ONLY: xp
#else
    use ondm_gen_com_m, ONLY: imm, A2cm, lperiod, bg, at, indi2
    use ondm_tab_imm_m, ONLY: iwmax2, xp
#endif
    use angular_functions
    use ml_in_ndm_module, ONLY: rangml, imm_neigh, &
                                weighted, desc_forces, linvisible
    use derived_types, only: config_real
    use module_neigh_local, only: r_cut 
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
    integer  :: iw, iw1, iw2

    integer  :: ia, ja, ia_n
    double precision     :: r2_ji, r_ji
    logical  :: desc_forces_local


    desc_forces_local = desc_forces .and. config_real(iconf)%has_force

    if ((i_start_at == 0) .and. (i_final_at == 0)) then
      d_n_neigh(:) = 0
      d_kind_neigh(:, :) = 0
      ! config_desc(iconf)%energy(:,:)=0.d0
      return
    end if

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
    d_n_neigh(:) = 0
    d_kind_neigh(:, :) = 0
    ! config_desc(iconf)%energy(:,:)=0.d0

#ifdef MLD_NDM
!!$    if (i_start_at == 1) iw2 = 0
!!$    if (i_start_at > 1) iw2 = iwmax2(i_start_at - 1)
    iw2=0
#else
    if (i_start_at == 1) iw2 = 0
    if (i_start_at > 1) iw2 = iwmax2(i_start_at - 1)
#endif
    do ja = i_start_at, i_final_at

      if (linvisible .and. weighted) then
        if (config_real(iconf)%invisible_per_type(config_real(iconf)%itype(ja))) cycle
      end if
      ! begin small box or not 1/
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

      ia_n = 0
      do iw = iw1, iw2
        !begin small box or not 2/
#ifdef MLD_NDM
!!$        if (small) then
#else
        if (small) then
#endif
          ia = config_real(iconf)%kind_neigh(ja, iw)
          !write (*,*) 'debug kind_neigh', ia
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
        !end small box or not 2/

        if (weighted) then
          if (linvisible) then
            if (config_real(iconf)%invisible_per_type(config_real(iconf)%itype(ia))) cycle
          end if
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

        if (r_ji .gt. r_cut) cycle
        ia_n = ia_n + 1
        d_kind_neigh(ja, ia_n) = ia
      end do                  ! iw
      d_n_neigh(ja) = ia_n
      ! read __desc energy__
      ! read __desc force__
      ! config_desc(iconf)%energy(1:dim_xdesc,ja) = local_desc(1:dim_xdesc,ja)
      if (desc_forces_local) then
        do ia = 0, ia_n
          ! config_desc(iconf)%force(1:dim_xdesc,ja, ia,1:3) = local_desc_deriv(1:dim_xdesc,ja,ia,1:3)
        end do
      end if
    end do                  ! ja
    deallocate (xpnp)

  end subroutine compute_dump_desc

end module
