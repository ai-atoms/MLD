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

! subroutine compute_mtp(i_start_at,i_final_at,d_n_neigh, d_kind_neigh, local_mtp_out,local_mtp_deriv_out, iconf)
module module_compute_mtp 
contains 
subroutine compute_mtp(i_start_at, i_final_at, d_n_neigh, d_kind_neigh, iconf)

  use module_kind_variables, ONLY: double
#ifdef MLD_NDM
  use gen_com_m, ONLY: A2cm, lperiod
  use gen_com_m_ml, ONLY: imm, bg, at
  use tab_imm_m_ml, ONLY: xp
#else
  use ondm_gen_com_m, ONLY: imm, A2cm, lperiod, bg, at, indi2
  use ondm_tab_imm_m, ONLY: iwmax2, xp
#endif
  use angular_functions
  use ml_in_ndm_module, ONLY: rangml, one_pi, imm_neigh, &
                              weighted, mtp_rad_order, &
                              mtp_poly_min, desc_forces, linvisible
  use module_neigh_local, only: r_cut 
  use derived_types, only: config_real, config_desc
  use time_check_general, only: time, tot_time, debug_time, MY_MPI_WTIME
#ifdef MLD_NDM
  use notperiod_mod
#else
  use ondm_transform_coord, only: ondm_notperiod
#endif

  implicit none

  integer, intent(in)  :: i_start_at, i_final_at
  integer, dimension(imm), intent(out)   :: d_n_neigh
  integer, dimension(imm, imm_neigh), intent(out)    :: d_kind_neigh
  !real(kind(0.d0)),dimension(mtp_dim,imm),intent(out) :: local_mtp_out
  !real(kind(0.d0)),dimension(mtp_dim,imm, 0:imm_neigh, 3),intent(out) :: local_mtp_deriv_out
  integer, optional    :: iconf
  logical  :: small

  real(kind(0.d0)), dimension(mtp_rad_order)   :: tmp_mtp1
  real(kind(0.d0)), dimension(mtp_rad_order, 3)      :: tmp_mtp2
  real(kind(0.d0)), dimension(mtp_rad_order, 3, 3)   :: tmp_mtp3
  !real(kind(0.d0)), dimension(mtp_rad_order, 3, 3, 3)      :: tmp_mtp4
  real(kind(0.d0)), dimension(mtp_rad_order, 27)      :: tmp_mtp4
  !real(kind(0.d0)), dimension(mtp_rad_order, 3, 3, 3, 3)   :: tmp_mtp5
  real(kind(0.d0)), dimension(mtp_rad_order, 81)   :: tmp_mtp5
  !real(kind(0.d0)), dimension(mtp_rad_order, 3, 3, 3, 3, 3)      :: tmp_mtp6
  real(kind(0.d0)), dimension(mtp_rad_order, 243)      :: tmp_mtp6
  real(kind(0.d0)), dimension(:, :, :), allocatable  :: d_tmp_mtp1
  real(kind(0.d0)), dimension(:, :, :, :), allocatable     :: d_tmp_mtp2
  real(kind(0.d0)), dimension(:, :, :, :, :), allocatable  :: d_tmp_mtp3
  !real(kind(0.d0)),dimension(:, :, :, :, :, :), allocatable     :: d_tmp_mtp4
  real(kind(0.d0)), dimension(:, :, :, : ), allocatable     :: d_tmp_mtp4
  !real(kind(0.d0)), dimension(:, :, :, :, :, :, :), allocatable  :: d_tmp_mtp5
  real(kind(0.d0)), dimension(:, :, :, :), allocatable  :: d_tmp_mtp5
  !real(kind(0.d0)), dimension(:, :, :, :, :, :, :, :), allocatable     :: d_tmp_mtp6
  real(kind(0.d0)), dimension(:, :, :, :), allocatable     :: d_tmp_mtp6

  real(kind(0.d0)), dimension(:, :), allocatable     :: vx_ji
  real(kind(0.d0)), dimension(:), allocatable  :: vr_ji


  real(kind(0.d0)), dimension(:, :), allocatable     :: tmpx_ji
  real(kind(0.d0)), dimension(:), allocatable  :: tmpr_ji


  real(kind(0.d0)), dimension(:, :, :), allocatable  :: f_vr_ji
  real(kind(0.d0)), dimension(:, :, :), allocatable  :: d_f_vr_ji

  real(double), dimension(:, :), allocatable   :: xpnp
  real(double), dimension(:), allocatable   :: ttmmpp, v_fcut, v_dfcut
  real(double), dimension(3) :: dxp_ji, ds
  real(kind(0.d0))     :: tmp1, tmp2, tmp3, tmp4, tmp5
  integer  :: iw, iw1, iw2

  integer  :: ip, ipo, ip1, ip2, ip3, ix
  integer  :: b1, b2, b3, b4, b5, nu, ibb3, ibb4, ibb5
  integer  :: ia, ja, ia_n, icount_e, icount_f
  double precision     :: r2_ji, r_ji, a_ip1_ip2, b_ip1_ip2, c_ip1_ip2
  double precision     :: factor_ia, factor_ja, fcut, dfcut
  ! double precision ::  tmp_local, tmp_local_3d(3)
  logical  :: desc_forces_local


  desc_forces_local = desc_forces .and. (config_real(iconf)%has_force .or. config_real(iconf)%has_stress)
  if (allocated(v_fcut))   deallocate(v_fcut) ;  allocate(v_fcut(imm_neigh))
  if (allocated(v_dfcut)) deallocate(v_dfcut) ; allocate(v_dfcut(imm_neigh))
  if (desc_forces_local) then
    if (allocated(d_tmp_mtp1)) deallocate (d_tmp_mtp1); allocate (d_tmp_mtp1(mtp_rad_order, imm_neigh, 3))
    if (allocated(d_tmp_mtp2)) deallocate (d_tmp_mtp2); allocate (d_tmp_mtp2(mtp_rad_order, 3, imm_neigh, 3))
    if (allocated(d_tmp_mtp3)) deallocate (d_tmp_mtp3); allocate (d_tmp_mtp3(mtp_rad_order, 3, 3, imm_neigh, 3))
    !if (allocated(d_tmp_mtp4)) deallocate (d_tmp_mtp4); allocate (d_tmp_mtp4(mtp_rad_order, 3, 3, 3, imm_neigh, 3))
    if (allocated(d_tmp_mtp4)) deallocate (d_tmp_mtp4); allocate (d_tmp_mtp4(mtp_rad_order, 27, imm_neigh, 3))
    !if (allocated(d_tmp_mtp5)) deallocate (d_tmp_mtp5); allocate (d_tmp_mtp5(mtp_rad_order, 3, 3, 3, 3, imm_neigh, 3))
    if (allocated(d_tmp_mtp5)) deallocate (d_tmp_mtp5); allocate (d_tmp_mtp5(mtp_rad_order, 81, imm_neigh, 3))
    !if (allocated(d_tmp_mtp6)) deallocate (d_tmp_mtp6); allocate (d_tmp_mtp6(mtp_rad_order, 3, 3, 3, 3, 3, imm_neigh, 3))
    if (allocated(d_tmp_mtp6)) deallocate (d_tmp_mtp6); allocate (d_tmp_mtp6(mtp_rad_order, 243, imm_neigh, 3))
  end if

  if ((i_start_at == 0) .and. (i_final_at == 0)) then
    d_n_neigh(:) = 0
    d_kind_neigh(:, :) = 0
    config_desc(iconf)%energy(:, :) = 0.d0
    ! local_mtp_deriv_out(:,:,:,:)=0.d0
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
  config_desc(iconf)%energy(:, :) = 0.d0
  if (desc_forces_local) config_desc(iconf)%force(:, :, :, :) = 0.d0

#ifdef MLD_NDM
  iw2=0
!!$  if (i_start_at == 1) iw2 = 0
!!$  if (i_start_at > 1) iw2 = iwmax2(i_start_at - 1)
#else
  if (i_start_at == 1) iw2 = 0
  if (i_start_at > 1) iw2 = iwmax2(i_start_at - 1)
#endif

  if (debug_time) time(1) = MY_MPI_WTIME()
  do ja = i_start_at, i_final_at

    if (debug_time) time(3) = MY_MPI_WTIME()
    if (allocated(tmpx_ji)) deallocate (tmpx_ji); allocate (tmpx_ji(3, imm_neigh))
    if (allocated(tmpr_ji)) deallocate (tmpr_ji); allocate (tmpr_ji(imm_neigh))

    ! begin small box or not 1/
    if (linvisible .and. weighted) then
      if (config_real(iconf)%invisible_per_type(config_real(iconf)%itype(ja))) cycle
    end if
    if (weighted) then
      factor_ja = config_real(iconf)%weight_per_type(config_real(iconf)%itype(ja))                       ! /factor_weight_mass
    else
      factor_ja = 1.d0
    end if
#ifdef MLD_NDM
!!$    if (small) then
#else
    if (small) then
#endif
      iw1 = 1
      iw2 = config_real(iconf)%n_neigh(ja)
      ! debug write (*,'("rangml, debug neigh",4i5)') rangml, ja, iw1,  iw2
      ! debug write (*,*) rangml, config_real(iconf)%filename
#ifdef MLD_NDM
!!$    else
!!$      iw1 = iw2 + 1
!!$      iw2 = iwmax2(ja)
!!$      ! write (*,*) 'debug neigh', iw1, iw2
!!$    end if
#else
    else
      iw1 = iw2 + 1
      iw2 = iwmax2(ja)
      ! write (*,*) 'debug neigh', iw1, iw2
    end if
#endif
    ! end small box or not 1/

    tmp_mtp1(:) = 0.d0
    if (desc_forces_local) d_tmp_mtp1(:, :, :) = 0.d0
    tmp_mtp2(:, :) = 0.d0
    if (desc_forces_local) d_tmp_mtp2(:, :, :, :) = 0.d0
    tmp_mtp3(:, :, :) = 0.d0
    if (desc_forces_local) d_tmp_mtp3(:, :, :, :, :) = 0.d0
    !tmp_mtp4(:, :, :, :) = 0.d0
    !if (desc_forces_local) d_tmp_mtp4(:, :, :, :, :, :) = 0.d0
    tmp_mtp4(:, :) = 0.d0
    if (desc_forces_local) d_tmp_mtp4(:, :, :, :) = 0.d0
    !tmp_mtp5(:, :, :, :, :) = 0.d0
    !if (desc_forces_local) d_tmp_mtp5(:, :, :, :, :, :, :) = 0.d0
    tmp_mtp5(:, :) = 0.d0
    if (desc_forces_local) d_tmp_mtp5(:, :, :, :) = 0.d0
    !tmp_mtp6(:, :, :, :, :, :) = 0.d0
    !if (desc_forces_local) d_tmp_mtp6(:, :, :, :, :, :, :, :) = 0.d0
    tmp_mtp6(:, :) = 0.d0
    if (desc_forces_local) d_tmp_mtp6(:, :, :, :) = 0.d0

    if (weighted) then
      factor_ja = config_real(iconf)%weight_per_type(config_real(iconf)%itype(ja))                       ! /factor_weight_mass
    else
      factor_ja = 1.d0
    end if

    if (debug_time) then
      time(4) = MY_MPI_WTIME()
      tot_time(2) = tot_time(2) + time(4) - time(3)
    end if

    ia_n = 0
    do iw = iw1, iw2
      ! begin small box or not 2/
#ifdef MLD_NDM
!!$      if (small) then
#else
      if (small) then
#endif
        ia = config_real(iconf)%kind_neigh(ja, iw)
        ! write (*,*) 'debug kind_neigh', ia
#ifdef MLD_NDM
!!$      else
!!$        ia = indi2(iw)
!!$        if (ja == ia) cycle
!!$      end if
#else
      else
        ia = indi2(iw)
        if (ja == ia) cycle
      end if
#endif
      ! end small box or not 2/

      if (linvisible .and. weighted) then
        if (config_real(iconf)%invisible_per_type(config_real(iconf)%itype(ia))) cycle
      end if
      if (weighted) then
        factor_ia = config_real(iconf)%weight_per_type(config_real(iconf)%itype(ia))                       ! /factor_weight_mass
      else
        factor_ia = 1.d0
      end if

#ifdef MLD_NDM
!!$      if (small) then
#else
      if (small) then
#endif
        r_ji = config_real(iconf)%r_ij(ja, iw)
        dxp_ji(:) = config_real(iconf)%u_ij(ja, iw, :)
#ifdef MLD_NDM
!!$      else
!!$        dxp_ji(1:3) = xpnp(1:3, ia) - xpnp(1:3, ja)
!!$        ds = MatMul(dxp_ji, bg)
!!$        WHERE ((ds .GT. 0.5d0) .OR. (ds .LT. -0.5d0))
!!$          ds(1:3) = ds(1:3) - Dble(Nint(ds(1:3)))
!!$        END WHERE
!!$        dxp_ji = MatMul(at, ds)/A2cm
!!$        r2_ji = Sum(dxp_ji(1:3)**2)
!!$        r_ji = dsqrt(r2_ji)
!!$      end if
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
      tmpx_ji(:, ia_n) = dxp_ji(:)
      tmpr_ji(ia_n) = r_ji
      fcut = 0.5d0*(cos(one_pi*r_ji/r_cut) + 1d0)*factor_ia
      dfcut = -0.5d0*one_pi/r_cut*sin(one_pi*r_ji/r_cut)*factor_ia
      v_fcut(ia_n) = fcut
      v_dfcut(ia_n) = dfcut/r_ji
      !v_fcut(ia_n) = 1.d0
      !v_dfcut(ia_n) = 0.d0
      !r_fcut(1:max_neigh_local) = ((r_central(1:max_neigh_local)/r_cut)**2 - 1)**2
      !if (desc_forces_local) then
      ! do ii = 1, 3
      !   d_r_fcut(ii, 1:max_neigh_local) = tmpcos_dxp(ii, 1:max_neigh_local)* &
      !                                  4.d0*((r_central(1:max_neigh_local)/r_cut)**2 - 1)*r_central(1:max_neigh_local)/r_cut**2
      !end do
      d_kind_neigh(ja, ia_n) = ia
    end do                  ! iw

    if (debug_time) then
      time(5) = MY_MPI_WTIME()
      tot_time(3) = tot_time(3) + time(5) - time(4)
    end if
    ! ia_n is the max of neighbours
    d_n_neigh(ja) = ia_n
    if (ia_n > imm_neigh) then
      write (6, *) 'Fatal for atom ja for which the number of neighbours is higher than the admitted MAXIMUM', ja, imm_neigh
      stop 'imm_neigh not large enough or the structure is wrong - atoms too close'
    end if

    if (allocated(vx_ji)) deallocate (vx_ji); allocate (vx_ji(3, ia_n))
    if (allocated(vr_ji)) deallocate (vr_ji); allocate (vr_ji(ia_n))
    vx_ji(1:3, 1:ia_n) = tmpx_ji(1:3, 1:ia_n)
    deallocate (tmpx_ji)
    vr_ji(1:ia_n) = tmpr_ji(1:ia_n)
    deallocate (tmpr_ji)
    if (allocated(f_vr_ji)) deallocate (f_vr_ji); allocate (f_vr_ji(mtp_rad_order, 0:5, ia_n))
    if (allocated(d_f_vr_ji)) deallocate (d_f_vr_ji); allocate (d_f_vr_ji(mtp_rad_order, 0:5, ia_n))
    if (allocated(ttmmpp)) deallocate(ttmmpp) ; allocate (ttmmpp(ia_n))
    ! initialize the function
    do ip = 1, mtp_rad_order
      ipo = -ip - mtp_poly_min + 1
      do nu = 0, 5
        ttmmpp(1:ia_n)  = vr_ji(1:ia_n)**(ipo - nu)
        f_vr_ji(ip, nu, 1:ia_n) =ttmmpp(1:ia_n)*v_fcut(1:ia_n)
        d_f_vr_ji(ip, nu, 1:ia_n) = dble(ipo-nu)*vr_ji(1:ia_n)**(ipo - 2 - nu)*v_fcut(1:ia_n) &
                                    + v_dfcut(1:ia_n)*ttmmpp(1:ia_n)
      end do
    end do

    do ip = 1, mtp_rad_order
      ipo = -ip - mtp_poly_min + 1

      ! zero order mtp(ipo, 0) ~ mtp(ipo)
      tmp_mtp1(ip) = SUM(f_vr_ji(ip, 0, :))
      if (desc_forces_local) then
        do ia = 1, ia_n
          d_tmp_mtp1(ip, ia, :) = d_f_vr_ji(ip, 0, ia)*vx_ji(:, ia)
        end do
      end if

      ! second order mtp(ip, 1) ~ mtp(ip, b1)
      do b1 = 1, 3
        do ia = 1, ia_n
          tmp_mtp2(ip, b1) = tmp_mtp2(ip, b1) + f_vr_ji(ip, 1, ia)*vx_ji(b1, ia)
        end do
      end do

      if (desc_forces_local) then
        do b1 = 1, 3
          do ia = 1, ia_n
            do ix = 1, 3
              tmp1 = 0.d0
              if (ix == b1) tmp1 = f_vr_ji(ip, 1, ia)
              d_tmp_mtp2(ip, b1, ia, ix) = d_f_vr_ji(ip, 1, ia)*vx_ji(ix, ia)*vx_ji(b1, ia) + tmp1
            end do
          end do
        end do
      end if

      if (debug_time) then
        time(6) = MY_MPI_WTIME()
        tot_time(4) = tot_time(4) + time(6) - time(5)
       end if
      ! third order mtp(ip, 2) ~ mtp(ip, b1, b2)
      do b1 = 1, 3
        do b2 = 1, 3
          do ia = 1, ia_n
            tmp_mtp3(ip, b1, b2) = tmp_mtp3(ip, b1, b2) + f_vr_ji(ip, 2, ia)*vx_ji(b1, ia)*vx_ji(b2, ia)
          end do
        end do
      end do

      if (desc_forces_local) then
        do b1 = 1, 3
          do b2 = 1, 3
            do ia = 1, ia_n
              do ix = 1, 3
                tmp1 = 0.d0
                tmp2 = 0.d0
                if (b1 == ix) tmp1 = f_vr_ji(ip, 2, ia)*vx_ji(b2, ia)
                if (b2 == ix) tmp2 = f_vr_ji(ip, 2, ia)*vx_ji(b1, ia)
                d_tmp_mtp3(ip, b1, b2, ia, ix) = d_f_vr_ji(ip, 2, ia)*vx_ji(ix, ia)*vx_ji(b1, ia)*vx_ji(b2, ia) + tmp1 + tmp2
              end do
            end do
          end do
        end do
      end if

      if (debug_time) then
        time(7) = MY_MPI_WTIME()
        tot_time(5) = tot_time(5) + time(7) - time(6)
       end if

      ! fourth order mtp(ip, 3) ~ mtp(ip, b1, b2,b3)
      do b1 = 1, 3
        do b2 = 1, 3
          do b3 = 1, 3
            ibb3 = (b1-1)*9 + (b2-1)*3 + b3
            do ia = 1, ia_n
              !tmp_mtp4(ip, b1, b2, b3) = tmp_mtp4(ip, b1, b2, b3) + f_vr_ji(ip, 3, ia)*vx_ji(b1, ia)*vx_ji(b2, ia)*vx_ji(b3, ia)
              tmp_mtp4(ip, ibb3) = tmp_mtp4(ip, ibb3) + f_vr_ji(ip, 3, ia)*vx_ji(b1, ia)*vx_ji(b2, ia)*vx_ji(b3, ia)
            end do
          end do
        end do
      end do

      if (desc_forces_local) then
        do b1 = 1, 3
          do b2 = 1, 3
            do b3 = 1, 3
              ibb3 = (b1-1)*9 + (b2-1)*3 + b3
              do ia = 1, ia_n
                do ix = 1, 3
                  tmp1 = 0.d0
                  tmp2 = 0.d0
                  tmp3 = 0.d0
                  if (b1 == ix) tmp1 = f_vr_ji(ip, 3, ia)*vx_ji(b2, ia)*vx_ji(b3, ia)
                  if (b2 == ix) tmp2 = f_vr_ji(ip, 3, ia)*vx_ji(b1, ia)*vx_ji(b3, ia)
                  if (b3 == ix) tmp3 = f_vr_ji(ip, 3, ia)*vx_ji(b1, ia)*vx_ji(b2, ia)
                  d_tmp_mtp4(ip, ibb3,  ia, ix) = d_f_vr_ji(ip, 3, ia)*vx_ji(ix, ia)*vx_ji(b1, ia)*vx_ji(b2, ia)*vx_ji(b3, ia) + tmp1 + tmp2 + tmp3
                end do
              end do
            end do
          end do
        end do
      end if

      if (debug_time) then
        time(8) = MY_MPI_WTIME()
        tot_time(6) = tot_time(6) + time(8) - time(7)
       end if


      ! fifth order
      ! mtp(ip, 4) ~ mtp(ip, b1, b2, b3, b4)
      ibb4 = 0
      do b1 = 1, 3
        do b2 = 1, 3
          do b3 = 1, 3
            do b4 = 1, 3
              ibb4 = ibb4 + 1
              do ia = 1, ia_n
                !tmp_mtp5(ip, b1, b2, b3, b4) = tmp_mtp5(ip, b1, b2, b3, b4) + f_vr_ji(ip, 4, ia)*vx_ji(b1, ia)*vx_ji(b2, ia)*vx_ji(b3, ia)*vx_ji(b4, ia)
                tmp_mtp5(ip, ibb4) = tmp_mtp5(ip, ibb4) + f_vr_ji(ip, 4, ia)*vx_ji(b1, ia)*vx_ji(b2, ia)*vx_ji(b3, ia)*vx_ji(b4, ia)
              end do
            end do
          end do
        end do
      end do

      if (desc_forces_local) then
        ibb4 = 0
        do b1 = 1, 3
          do b2 = 1, 3
            do b3 = 1, 3
              do b4 = 1, 3
                ibb4 = ibb4 + 1
                do ia = 1, ia_n
                  do ix = 1, 3
                    tmp1 = 0.d0
                    tmp2 = 0.d0
                    tmp3 = 0.d0
                    tmp4 = 0.d0
                    !if (b1 == ix) tmp1 = f_vr_ji(ip, 4, ia)*vx_ji(b2, ia)*vx_ji(b3, ia)*vx_ji(b4, ia)
                    !if (b2 == ix) tmp2 = f_vr_ji(ip, 4, ia)*vx_ji(b1, ia)*vx_ji(b3, ia)*vx_ji(b4, ia)
                    !if (b3 == ix) tmp3 = f_vr_ji(ip, 4, ia)*vx_ji(b1, ia)*vx_ji(b2, ia)*vx_ji(b4, ia)
                    !if (b4 == ix) tmp4 = f_vr_ji(ip, 4, ia)*vx_ji(b1, ia)*vx_ji(b2, ia)*vx_ji(b3, ia)
                    !d_tmp_mtp5(ip, b1, b2, b3, b4, ia, ix) = d_f_vr_ji(ip, 4, ia)*vx_ji(ix, ia)*vx_ji(b1, ia)*vx_ji(b2, ia)*vx_ji(b3, ia)*vx_ji(b4, ia) + tmp1 + tmp2 + tmp3 + tmp4
                    if (b1 == ix) tmp1 = f_vr_ji(ip, 4, ia)*vx_ji(b2, ia)*vx_ji(b3, ia)*vx_ji(b4, ia)
                    if (b2 == ix) tmp2 = f_vr_ji(ip, 4, ia)*vx_ji(b1, ia)*vx_ji(b3, ia)*vx_ji(b4, ia)
                    if (b3 == ix) tmp3 = f_vr_ji(ip, 4, ia)*vx_ji(b1, ia)*vx_ji(b2, ia)*vx_ji(b4, ia)
                    if (b4 == ix) tmp4 = f_vr_ji(ip, 4, ia)*vx_ji(b1, ia)*vx_ji(b2, ia)*vx_ji(b3, ia)
                    d_tmp_mtp5(ip, ibb4, ia, ix) = d_f_vr_ji(ip, 4, ia)*vx_ji(ix, ia)*vx_ji(b1, ia)*vx_ji(b2, ia)*vx_ji(b3, ia)*vx_ji(b4, ia) + tmp1 + tmp2 + tmp3 + tmp4
                  end do
                end do
              end do
            end do
          end do
        end do
      end if


      if (debug_time) then
        time(9) = MY_MPI_WTIME()
        tot_time(7) = tot_time(7) + time(9) - time(8)
       end if

      ! sixth order mtp(ip, 5) ~ mtp(ip, b1, b2, b3, b4, b5)
     ibb5 = 0
      do b1 = 1, 3
        do b2 = 1, 3
          do b3 = 1, 3
            do b4 = 1, 3
              do b5 = 1, 3
                ibb5 = ibb5 + 1
                do ia = 1, ia_n
                  tmp_mtp6(ip, ibb5) = tmp_mtp6(ip, ibb5) + f_vr_ji(ip, 5, ia)*vx_ji(b1, ia)*vx_ji(b2, ia)*vx_ji(b3, ia)*vx_ji(b4, ia)*vx_ji(b5, ia)
                end do
              end do
            end do
          end do
        end do
      end do

      if (desc_forces_local) then
       ibb5 = 0
        do b1 = 1, 3
          do b2 = 1, 3
            do b3 = 1, 3
              do b4 = 1, 3
                do b5 = 1, 3
                  ibb5 = ibb5 + 1
                  do ia = 1, ia_n
                    do ix = 1, 3
                      tmp1 = 0.d0
                      tmp2 = 0.d0
                      tmp3 = 0.d0
                      tmp4 = 0.d0
                      tmp5 = 0.d0
                      if (b1 == ix) tmp1 = f_vr_ji(ip, 5, ia)*vx_ji(b2, ia)*vx_ji(b3, ia)*vx_ji(b4, ia)*vx_ji(b5, ia)
                      if (b2 == ix) tmp2 = f_vr_ji(ip, 5, ia)*vx_ji(b1, ia)*vx_ji(b3, ia)*vx_ji(b4, ia)*vx_ji(b5, ia)
                      if (b3 == ix) tmp3 = f_vr_ji(ip, 5, ia)*vx_ji(b1, ia)*vx_ji(b2, ia)*vx_ji(b4, ia)*vx_ji(b5, ia)
                      if (b4 == ix) tmp4 = f_vr_ji(ip, 5, ia)*vx_ji(b1, ia)*vx_ji(b2, ia)*vx_ji(b3, ia)*vx_ji(b5, ia)
                      if (b5 == ix) tmp5 = f_vr_ji(ip, 5, ia)*vx_ji(b1, ia)*vx_ji(b2, ia)*vx_ji(b3, ia)*vx_ji(b4, ia)
                      d_tmp_mtp6(ip, ibb5, ia, ix) = d_f_vr_ji(ip, 5, ia)*vx_ji(ix, ia)*vx_ji(b1, ia)*vx_ji(b2, ia)*vx_ji(b3, ia)*vx_ji(b4, ia)*vx_ji(b5, ia) &
                                                                   + tmp1 + tmp2 + tmp3 + tmp4 + tmp5
                    end do
                  end do
                end do
              end do
            end do
          end do
        end do
      end if

    end do                  ! first ip.

      if (debug_time) then
        time(10) = MY_MPI_WTIME()
        tot_time(8) = tot_time(8) + time(10) - time(9)
       end if
    ! compute the descriptor the reduction of the above tensors
    !01---------------------------------!
    ! alpha = ( \mu  )                  !
    !-----------------------------------!
    config_desc(iconf)%energy(1:mtp_rad_order, ja) = tmp_mtp1(1:mtp_rad_order)
    if (desc_forces_local) config_desc(iconf)%force(1:mtp_rad_order, ja, 1:ia_n, 1:3) = d_tmp_mtp1(1:mtp_rad_order, 1:ia_n, 1:3)

    !02---------------------------------!
    ! alpha = ( \mu_1      1)           !
    !         (     1  \mu_2)           !
    !             -----                 !
    ! yields \mu1*\mu2 components       !
    !-----------------------------------!
    config_desc(iconf)%energy(mtp_rad_order + 1:mtp_rad_order + mtp_rad_order**2, ja) = RESHAPE(MATMUL(tmp_mtp2, TRANSPOSE(tmp_mtp2)), (/mtp_rad_order**2/))
    icount_e = mtp_rad_order + mtp_rad_order**2 + 1

    if (desc_forces_local) then
      icount_f = mtp_rad_order + 1
      do ip1 = 1, mtp_rad_order
        do ip2 = 1, mtp_rad_order
          do b1 = 1, 3
            do ia = 1, ia_n
            do ix = 1, 3
              a_ip1_ip2 = d_tmp_mtp2(ip1, b1, ia, ix)*tmp_mtp2(ip2, b1)
              b_ip1_ip2 = tmp_mtp2(ip1, b1)*d_tmp_mtp2(ip2, b1, ia, ix)
              config_desc(iconf)%force(icount_f, ja, ia, ix) = config_desc(iconf)%force(icount_f, ja, ia, ix) + a_ip1_ip2 + b_ip1_ip2
              ! config_desc(iconf)%force(icount,ja,ia,iy) =  a_ip1_ip2  + b_ip1_ip2
            end do
            end do
          end do
          icount_f = icount_f + 1
        end do
      end do
    end if

    !02a--------------------------------!
    ! alpha = ( \mu_1      2)           !
    !         (     2  \mu_2)           !
    !             -----                 !
    ! yields \mu1*\mu2 components       !
    !-----------------------------------!
    do ip1 = 1, mtp_rad_order
      do ip2 = 1, mtp_rad_order
        do b1 = 1, 3
          do b2 = 1, 3
            ! config_desc(iconf)%energy(icount,ja)= config_desc(iconf)%energy(icount,ja) + tmp_mtp3(ip1,ix,iy)*tmp_mtp3(ip2,ix,iy)
            config_desc(iconf)%energy(icount_e, ja) = config_desc(iconf)%energy(icount_e, ja) + tmp_mtp3(ip1, b1, b2)*tmp_mtp3(ip2, b1, b2)
          end do
        end do
        icount_e = icount_e + 1
      end do
    end do

    if (desc_forces_local) then
      do ip1 = 1, mtp_rad_order
        do ip2 = 1, mtp_rad_order
          do b1 = 1, 3
            do b2 = 1, 3
              do ia = 1, ia_n
                do ix = 1, 3
                  a_ip1_ip2 = d_tmp_mtp3(ip1, b1, b2, ia, ix)*tmp_mtp3(ip2, b1, b2)
                  b_ip1_ip2 = tmp_mtp3(ip1, b1, b2)*d_tmp_mtp3(ip2, b1, b2, ia, ix)
                  config_desc(iconf)%force(icount_f, ja, ia, ix) = config_desc(iconf)%force(icount_f, ja, ia, ix) + a_ip1_ip2 + b_ip1_ip2
                  ! config_desc(iconf)%force(icount,ja,ia,iz) =  a_ip1_ip2 + b_ip1_ip2
                end do
              end do
            end do
          end do
          icount_f = icount_f + 1
        end do
      end do
    end if


    !02b--------------------------------!
    ! alpha = ( \mu_1      3)           !
    !         (     3  \mu_2)           !
    !             -----                 !
    ! yields \mu1*\mu2 components       !
    !-----------------------------------!
    do ip1 = 1, mtp_rad_order
      do ip2 = 1, mtp_rad_order
        do b1 = 1, 3
          do b2 = 1, 3
            do b3 = 1, 3
              ibb3 = (b1-1)*9 + (b2-1)*3 + b3
              !config_desc(iconf)%energy(icount_e, ja) = config_desc(iconf)%energy(icount_e, ja) + tmp_mtp4(ip1, b1, b2, b3)*tmp_mtp4(ip2, b1, b2, b3)
              config_desc(iconf)%energy(icount_e, ja) = config_desc(iconf)%energy(icount_e, ja) + tmp_mtp4(ip1, ibb3)*tmp_mtp4(ip2, ibb3)
            end do
          end do
        end do
        icount_e = icount_e + 1
      end do
    end do

    if (desc_forces_local) then
      do ip1 = 1, mtp_rad_order
        do ip2 = 1, mtp_rad_order
          do b1 = 1, 3
            do b2 = 1, 3
              do b3 = 1, 3
                ibb3 = (b1-1)*9 + (b2-1)*3 + b3
                do ia = 1, ia_n
                  do ix = 1, 3
                    a_ip1_ip2 = d_tmp_mtp4(ip1, ibb3, ia, ix)*tmp_mtp4(ip2, ibb3)
                    b_ip1_ip2 = tmp_mtp4(ip1, ibb3)*d_tmp_mtp4(ip2, ibb3, ia, ix)
                    config_desc(iconf)%force(icount_f, ja, ia, ix) = config_desc(iconf)%force(icount_f, ja, ia, ix) + a_ip1_ip2 + b_ip1_ip2
                    ! config_desc(iconf)%force(icount,ja,ia,iz) =  a_ip1_ip2 + b_ip1_ip2
                    !a_ip1_ip2 = d_tmp_mtp4(ip1, b1, b2, b3, ia, ix)*tmp_mtp4(ip2, b1, b2, b3)
                    !b_ip1_ip2 = tmp_mtp4(ip1, b1, b2, b3)*d_tmp_mtp4(ip2, b1, b2, b3, ia, ix)
                    !config_desc(iconf)%force(icount_f, ja, ia, ix) = config_desc(iconf)%force(icount_f, ja, ia, ix) + a_ip1_ip2 + b_ip1_ip2
                    !! config_desc(iconf)%force(icount,ja,ia,iz) =  a_ip1_ip2 + b_ip1_ip2
                  end do
                end do
              end do
            end do
          end do
          icount_f = icount_f + 1
        end do
      end do
    end if


    !02c--------------------------------!
    ! alpha = ( \mu_1      4)           !
    !         (     4  \mu_2)           !
    !             -----                 !
    ! yields \mu1*\mu2 components       !
    !-----------------------------------!
    do ip1 = 1, mtp_rad_order
      do ip2 = 1, mtp_rad_order
        ibb4 = 0
        do b1 = 1, 3
          do b2 = 1, 3
            do b3 = 1, 3
              do b4 = 1, 3
                ibb4 = ibb4 + 1
                !config_desc(iconf)%energy(icount_e, ja) = config_desc(iconf)%energy(icount_e, ja) + tmp_mtp5(ip1, b1, b2, b3, b4)*tmp_mtp5(ip2, b1, b2, b3, b4)
                config_desc(iconf)%energy(icount_e, ja) = config_desc(iconf)%energy(icount_e, ja) + tmp_mtp5(ip1, ibb4)*tmp_mtp5(ip2, ibb4)
              end do
            end do
          end do
        end do
        icount_e = icount_e + 1
      end do
    end do

    if (desc_forces_local) then
      do ip1 = 1, mtp_rad_order
        do ip2 = 1, mtp_rad_order
          ibb4 = 0
          do b1 = 1, 3
            do b2 = 1, 3
              do b3 = 1, 3
                do b4 = 1, 3
                  ibb4 = ibb4 + 1
                  do ia = 1, ia_n
                    do ix = 1, 3
                      !a_ip1_ip2 = d_tmp_mtp5(ip1, b1, b2, b3, b4, ia, ix)*tmp_mtp5(ip2, b1, b2, b3, b4)
                      !b_ip1_ip2 = tmp_mtp5(ip1, b1, b2, b3, b4)*d_tmp_mtp5(ip2, b1, b2, b3, b4, ia, ix)
                      !config_desc(iconf)%force(icount_f, ja, ia, ix) = config_desc(iconf)%force(icount_f, ja, ia, ix) + a_ip1_ip2 + b_ip1_ip2

                      a_ip1_ip2 = d_tmp_mtp5(ip1, ibb4, ia, ix)*tmp_mtp5(ip2, ibb4)
                      b_ip1_ip2 = tmp_mtp5(ip1, ibb4)*d_tmp_mtp5(ip2, ibb4, ia, ix)
                      config_desc(iconf)%force(icount_f, ja, ia, ix) = config_desc(iconf)%force(icount_f, ja, ia, ix) + a_ip1_ip2 + b_ip1_ip2
                    end do
                  end do
                end do
              end do
            end do
          end do
          icount_f = icount_f + 1
        end do
      end do
    end if


    !02d--------------------------------!
    ! alpha = ( \mu_1      5)           !
    !         (     5  \mu_2)           !
    !             -----                 !
    ! yields \mu1*\mu2 components       !
    !-----------------------------------!
    do ip1 = 1, mtp_rad_order
      do ip2 = 1, mtp_rad_order
        ibb5 = 0
        do b1 = 1, 3
          do b2 = 1, 3
            do b3 = 1, 3
              do b4 = 1, 3
                do b5 = 1, 3
                  ibb5 = ibb5 + 1
                  !config_desc(iconf)%energy(icount_e, ja) = config_desc(iconf)%energy(icount_e, ja) + tmp_mtp6(ip1, b1, b2, b3, b4, b5)*tmp_mtp6(ip2, b1, b2, b3, b4, b5)
                  config_desc(iconf)%energy(icount_e, ja) = config_desc(iconf)%energy(icount_e, ja) + tmp_mtp6(ip1, ibb5)*tmp_mtp6(ip2, ibb5)
                end do
              end do
            end do
          end do
        end do
        icount_e = icount_e + 1
      end do
    end do

    if (desc_forces_local) then
      do ip1 = 1, mtp_rad_order
        do ip2 = 1, mtp_rad_order
          ibb5 = 0
          do b1 = 1, 3
            do b2 = 1, 3
              do b3 = 1, 3
                do b4 = 1, 3
                  do b5 = 1, 3
                    ibb5 = ibb5 + 1
                    do ia = 1, ia_n
                      do ix = 1, 3
                        !a_ip1_ip2 = d_tmp_mtp6(ip1, b1, b2, b3, b4, b5, ia, ix)*tmp_mtp6(ip2, b1, b2, b3, b4, b5)
                        !b_ip1_ip2 = tmp_mtp6(ip1, b1, b2, b3, b4, b5)*d_tmp_mtp6(ip2, b1, b2, b3, b4, b5, ia, ix)
                        !config_desc(iconf)%force(icount_f, ja, ia, ix) = config_desc(iconf)%force(icount_f, ja, ia, ix) + a_ip1_ip2 + b_ip1_ip2
                        a_ip1_ip2 = d_tmp_mtp6(ip1, ibb5, ia, ix)*tmp_mtp6(ip2, ibb5)
                        b_ip1_ip2 = tmp_mtp6(ip1, ibb5)*d_tmp_mtp6(ip2, ibb5, ia, ix)
                        config_desc(iconf)%force(icount_f, ja, ia, ix) = config_desc(iconf)%force(icount_f, ja, ia, ix) + a_ip1_ip2 + b_ip1_ip2
                      end do
                    end do
                  end do
                end do
              end do
            end do
          end do
          icount_f = icount_f + 1
        end do
      end do
    end if

    !03---------------------------------!
    ! alpha = ( \mu_1      1      1)    !
    !         (     1  \mu_2      0)    !
    !         (     1      0  \mu_3)    !
    !             -----                 !
    ! yields \mu1*\mu2*\mu3 components  !
    !-----------------------------------!
    do ip1 = 1, mtp_rad_order
      do ip2 = 1, mtp_rad_order
        do ip3 = 1, mtp_rad_order
          do b1 = 1, 3
            do b2 = 1, 3
              config_desc(iconf)%energy(icount_e, ja) = config_desc(iconf)%energy(icount_e, ja) + tmp_mtp2(ip1, b1)*tmp_mtp3(ip2, b1, b2)*tmp_mtp2(ip3, b2)
            end do
          end do
          icount_e = icount_e + 1
        end do
      end do
    end do


    if (desc_forces_local) then
      do ip1 = 1, mtp_rad_order
        do ip2 = 1, mtp_rad_order
          do ip3 = 1, mtp_rad_order
            do b1 = 1, 3
              do b2 = 1, 3
                do ia = 1, ia_n
                  do ix = 1, 3
                    a_ip1_ip2 = tmp_mtp2(ip1, b1)*d_tmp_mtp3(ip2, b1, b2, ia, ix)*tmp_mtp2(ip3, b2)
                    b_ip1_ip2 = d_tmp_mtp2(ip1, b1, ia, ix)*tmp_mtp3(ip2, b1, b2)*tmp_mtp2(ip3, b2)
                    c_ip1_ip2 = tmp_mtp2(ip1, b1)*tmp_mtp3(ip2, b1, b2)*d_tmp_mtp2(ip3, b2, ia, ix)
                    config_desc(iconf)%force(icount_f, ja, ia, ix) = config_desc(iconf)%force(icount_f, ja, ia, ix) + a_ip1_ip2 + b_ip1_ip2 + c_ip1_ip2
                    ! config_desc(iconf)%force(icount,ja,ia,iz) =  a_ip1_ip2 + b_ip1_ip2 + c_ip1_ip2
                  end do
                end do
              end do
            end do
            icount_f = icount_f + 1
          end do
        end do
      end do
    end if


    !03a--------------------------------!
    ! alpha = ( \mu_1      1      2)    !
    !         (     1  \mu_2      0)    !
    !         (     2      0  \mu_3)    !
    !             -----                 !
    ! yields \mu1*\mu2*\mu3 components  !
    !-----------------------------------!
    do ip1 = 1, mtp_rad_order
      do ip2 = 1, mtp_rad_order
        do ip3 = 1, mtp_rad_order
          do b1 = 1, 3
            do b2 = 1, 3
              do b3 = 1, 3
                ibb3 = (b1-1)*9 + (b2-1)*3 + b3
                config_desc(iconf)%energy(icount_e, ja) = config_desc(iconf)%energy(icount_e, ja) + tmp_mtp2(ip1, b1)*tmp_mtp4(ip2, ibb3)*tmp_mtp3(ip3, b2, b3)
                !config_desc(iconf)%energy(icount_e, ja) = config_desc(iconf)%energy(icount_e, ja) + tmp_mtp2(ip1, b1)*tmp_mtp4(ip2, b1, b2, b3)*tmp_mtp3(ip3, b2, b3)
              end do
            end do
          end do
          icount_e = icount_e + 1
        end do
      end do
    end do


    if (desc_forces_local) then
      do ip1 = 1, mtp_rad_order
        do ip2 = 1, mtp_rad_order
          do ip3 = 1, mtp_rad_order
            ibb3 = 0
            do b1 = 1, 3
              do b2 = 1, 3
                do b3 = 1, 3
                  ibb3 = ibb3 + 1
                  do ia = 1, ia_n
                    do ix = 1, 3
                      a_ip1_ip2 = tmp_mtp2(ip1, b1)*d_tmp_mtp4(ip2, ibb3,  ia, ix)*tmp_mtp3(ip3, b2, b3)
                      b_ip1_ip2 = d_tmp_mtp2(ip1, b1, ia, ix)*tmp_mtp4(ip2, ibb3)*tmp_mtp3(ip3, b2, b3)
                      c_ip1_ip2 = tmp_mtp2(ip1, b1)*tmp_mtp4(ip2, ibb3)*d_tmp_mtp3(ip3, b2, b3, ia, ix)
                      config_desc(iconf)%force(icount_f, ja, ia, ix) = config_desc(iconf)%force(icount_f, ja, ia, ix) + a_ip1_ip2 + b_ip1_ip2 + c_ip1_ip2
                      ! config_desc(iconf)%force(icount,ja,ia,iz) =  a_ip1_ip2 + b_ip1_ip2 + c_ip1_ip2
                      !a_ip1_ip2 = tmp_mtp2(ip1, b1)*d_tmp_mtp4(ip2, b1, b2, b3, ia, ix)*tmp_mtp3(ip3, b2, b3)
                      !b_ip1_ip2 = d_tmp_mtp2(ip1, b1, ia, ix)*tmp_mtp4(ip2, b1, b2, b3)*tmp_mtp3(ip3, b2, b3)
                      !c_ip1_ip2 = tmp_mtp2(ip1, b1)*tmp_mtp4(ip2, b1, b2, b3)*d_tmp_mtp3(ip3, b2, b3, ia, ix)
                      !config_desc(iconf)%force(icount_f, ja, ia, ix) = config_desc(iconf)%force(icount_f, ja, ia, ix) + a_ip1_ip2 + b_ip1_ip2 + c_ip1_ip2
                      !! config_desc(iconf)%force(icount,ja,ia,iz) =  a_ip1_ip2 + b_ip1_ip2 + c_ip1_ip2
                    end do
                  end do
                end do
              end do
            end do
            icount_f = icount_f + 1
          end do
        end do
      end do
    end if


    !03b--------------------------------!
    ! alpha = ( \mu_1      1      3)    !
    !         (     1  \mu_2      0)    !
    !         (     3      0  \mu_3)    !
    !             -----                 !
    ! yields \mu1*\mu2*\mu3 components  !
    !-----------------------------------!
    do ip1 = 1, mtp_rad_order
      do ip2 = 1, mtp_rad_order
        do ip3 = 1, mtp_rad_order
          ibb4 = 0
          do b1 = 1, 3
            ibb3 = 0
            do b2 = 1, 3
              do b3 = 1, 3
                do b4 = 1, 3
                  ibb4 = ibb4 + 1
                  ibb3 = ibb3 + 1
                  !ibb3 = (b2-1)*9 + (b3-1)*3 + b4
                  config_desc(iconf)%energy(icount_e, ja) = config_desc(iconf)%energy(icount_e, ja) + tmp_mtp5(ip1, ibb4)*tmp_mtp2(ip2, b1)*tmp_mtp4(ip3, ibb3)
                  !config_desc(iconf)%energy(icount_e, ja) = config_desc(iconf)%energy(icount_e, ja) + tmp_mtp5(ip1, b1, b2, b3, b4)*tmp_mtp2(ip2, b1)*tmp_mtp4(ip3, b2, b3, b4)
                end do
              end do
            end do
          end do
          icount_e = icount_e + 1
        end do
      end do
    end do


    if (desc_forces_local) then
      do ip1 = 1, mtp_rad_order
        do ip2 = 1, mtp_rad_order
          do ip3 = 1, mtp_rad_order
            ibb4 = 0
            do b1 = 1, 3
              ibb3 = 0
              do b2 = 1, 3
                do b3 = 1, 3
                  do b4 = 1, 3
                    !ibb3 =  (b2-1)*9  + (b3-1)*3 + b4
                    ibb3 = ibb3 + 1
                    ibb4 = ibb4 + 1
                    do ia = 1, ia_n
                      do ix = 1, 3
                        a_ip1_ip2 = tmp_mtp5(ip1, ibb4)*d_tmp_mtp2(ip2, b1, ia, ix)*tmp_mtp4(ip3, ibb3)
                        b_ip1_ip2 = d_tmp_mtp5(ip1, ibb4, ia, ix)*tmp_mtp2(ip2, b1)*tmp_mtp4(ip3, ibb3)
                        c_ip1_ip2 = tmp_mtp5(ip1, ibb4)*tmp_mtp2(ip2, b1)*d_tmp_mtp4(ip3, ibb3, ia, ix)
                        config_desc(iconf)%force(icount_f, ja, ia, ix) = config_desc(iconf)%force(icount_f, ja, ia, ix) + a_ip1_ip2 + b_ip1_ip2 + c_ip1_ip2
                        ! config_desc(iconf)%force(icount,ja,ia,iz) =  a_ip1_ip2 + b_ip1_ip2 + c_ip1_ip2
                        !a_ip1_ip2 = tmp_mtp5(ip1, b1, b2, b3, b4)*d_tmp_mtp2(ip2, b1, ia, ix)*tmp_mtp4(ip3, b2, b3, b4)
                        !b_ip1_ip2 = d_tmp_mtp5(ip1, b1, b2, b3, b4, ia, ix)*tmp_mtp2(ip2, b1)*tmp_mtp4(ip3, b2, b3, b4)
                        !c_ip1_ip2 = tmp_mtp5(ip1, b1, b2, b3, b4)*tmp_mtp2(ip2, b1)*d_tmp_mtp4(ip3, b2, b3, b4, ia, ix)
                        !config_desc(iconf)%force(icount_f, ja, ia, ix) = config_desc(iconf)%force(icount_f, ja, ia, ix) + a_ip1_ip2 + b_ip1_ip2 + c_ip1_ip2
                        !! config_desc(iconf)%force(icount,ja,ia,iz) =  a_ip1_ip2 + b_ip1_ip2 + c_ip1_ip2
                      end do
                    end do
                  end do
                end do
              end do
            end do
            icount_f = icount_f + 1
          end do
        end do
      end do
    end if


    !03c--------------------------------!
    ! alpha = ( \mu_1      1      4)    !
    !         (     1  \mu_2      0)    !
    !         (     4      0  \mu_3)    !
    !             -----                 !
    ! yields \mu1*\mu2*\mu3 components  !
    !-----------------------------------!
    do ip1 = 1, mtp_rad_order
      do ip2 = 1, mtp_rad_order
        do ip3 = 1, mtp_rad_order
          ibb5 = 0
          do b1 = 1, 3
            ibb4 = 0
            do b2 = 1, 3
              do b3 = 1, 3
                do b4 = 1, 3
                  do b5 = 1, 3
                    ibb5 = ibb5 + 1
                    ibb4 = ibb4 + 1
                    config_desc(iconf)%energy(icount_e, ja) = config_desc(iconf)%energy(icount_e, ja) + tmp_mtp6(ip1, ibb5)*tmp_mtp2(ip2, b1)*tmp_mtp5(ip3, ibb4)
                    !config_desc(iconf)%energy(icount_e, ja) = config_desc(iconf)%energy(icount_e, ja) + tmp_mtp6(ip1, b1, b2, b3, b4, b5)*tmp_mtp2(ip2, b1)*tmp_mtp5(ip3, b2, b3, b4, b5)
                  end do
                end do
              end do
            end do
          end do
          icount_e = icount_e + 1
        end do
      end do
    end do


    if (desc_forces_local) then
      do ip1 = 1, mtp_rad_order
        do ip2 = 1, mtp_rad_order
          do ip3 = 1, mtp_rad_order
            ibb5 = 0
            do b1 = 1, 3
              ibb4 = 0
              do b2 = 1, 3
                do b3 = 1, 3
                  do b4 = 1, 3
                    do b5 = 1, 3
                      ibb5 = ibb5 + 1
                      ibb4 = ibb4 + 1
                      do ia = 1, ia_n
                        do ix = 1, 3
                          !config_desc(iconf)%force(icount_f, ja, ia, ix) = config_desc(iconf)%force(icount_f, ja, ia, ix) + &
                          !                                                 tmp_mtp6(ip1, ibb5)*d_tmp_mtp2(ip2, b1, ia, ix)*tmp_mtp5(ip3, b2, b3, b4, b5) + &
                          !                                                 d_tmp_mtp6(ip1, b1, b2, b3, b4, b5, ia, ix)*tmp_mtp2(ip2, b1)*tmp_mtp5(ip3, b2, b3, b4, b5) + &
                          !                                                 tmp_mtp6(ip1, b1, b2, b3, b4, b5)*tmp_mtp2(ip2, b1)*d_tmp_mtp5(ip3, b2, b3, b4, b5, ia, ix)
                          config_desc(iconf)%force(icount_f, ja, ia, ix) = config_desc(iconf)%force(icount_f, ja, ia, ix) + &
                                                                           tmp_mtp6(ip1, ibb5)*d_tmp_mtp2(ip2, b1, ia, ix)*tmp_mtp5(ip3, ibb4) + &
                                                                           d_tmp_mtp6(ip1, ibb5, ia, ix)*tmp_mtp2(ip2, b1)*tmp_mtp5(ip3, ibb4) + &
                                                                           tmp_mtp6(ip1, ibb5)*tmp_mtp2(ip2, b1)*d_tmp_mtp5(ip3, ibb4, ia, ix)
                        end do
                      end do
                    end do
                  end do
                end do
              end do
            end do
            icount_f = icount_f + 1
          end do
        end do
      end do
    end if

    if (debug_time) then
      time(11) = MY_MPI_WTIME()
      tot_time(9) = tot_time(9) + time(11) - time(10)
    end if
    if (desc_forces_local) then
      do ia = 1, ia_n
        config_desc(iconf)%force(:, ja, 0, :) = config_desc(iconf)%force(:, ja, 0, :) - config_desc(iconf)%force(:, ja, ia, :)
      end do                  ! ia from ia_m
    end if
    ! NON NORMALIZED VERSION
    config_desc(iconf)%energy(:, ja) = config_desc(iconf)%energy(:, ja)*factor_ja                      ! /dble(ia_n)
    ! debug if (rangml==0) write (6,*) local_mtp_out(1,ja), ja
    if (desc_forces_local) config_desc(iconf)%force(:, ja, :, :) = config_desc(iconf)%force(:, ja, :, :)*factor_ja              ! /dble(ia_n)
  end do                  ! ja

  if (debug_time) then
    time(2) = MY_MPI_WTIME()
    tot_time(1) = tot_time(1) + time(2) - time(1)
  end if

  if (debug_time) then
    call repport_time(2, 0.d0, tot_time(1), "ML: full")
    call repport_time(2, 0.d0, tot_time(2), "ML: main loop: init before nn  ")
    call repport_time(2, 0.d0, tot_time(3), "ML: main loop: nn loop ")
    call repport_time(4, 0.d0, tot_time(4), "ML: tensor 1 and 2 ")
    call repport_time(4, 0.d0, tot_time(5), "ML: tensor 3 ")
    call repport_time(4, 0.d0, tot_time(6), "ML: tensor 4 ")
    call repport_time(4, 0.d0, tot_time(7), "ML: tensor 5 ")
    call repport_time(4, 0.d0, tot_time(8), "ML: tensor 6 ")
    call repport_time(4, 0.d0, tot_time(9), "ML: full reduction ")
  end if




  deallocate (xpnp)

end subroutine compute_mtp
end module module_compute_mtp


subroutine gen_dimension_for_mtp()
  use ml_in_ndm_module, ONLY: rangml, mtp_dim, mtp_poly_min, mtp_poly_max, mtp_rad_order
  implicit none

  if (mtp_poly_max <= mtp_poly_min) then
    if (rangml == 0) write (6, *) 'Error in setting MTP descriptor. mpt_poly_max sould be larger than mtp_poly_min'
    stop 'error in MPT in gen_dimension_for_mtp'
  end if
  mtp_rad_order = mtp_poly_max - mtp_poly_min + 1
  mtp_dim = mtp_rad_order + 5*mtp_rad_order**2 + 4*mtp_rad_order**3
  !6 mtp_dim=mtp_rad_order + 5*mtp_rad_order**2 + 2*mtp_rad_order**3
  if (rangml == 0) write (6, *) 'ML: dimension of the descriptor space: ', mtp_dim
end subroutine gen_dimension_for_mtp


!real(kind=kind(1.d0)) function mtp_radial_tensor(type_mtp, r, ipo) result(func, d_func)
!integer :: type_mtp,ipo
!real(kind=kind(1.d0)) :: r, func, d_func
!
!  select  case (type_mtp)
!
!  case(1)
!    func= r**ipo
!    d_func= dble(ipo)*r**(ipo-2)
!  end select
!
!end function mtp_radial_tensor
