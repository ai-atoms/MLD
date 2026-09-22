
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

! subroutine compute_pow_so4(i_start_at,i_final_at,d_n_neigh, d_kind_neigh, local_pow_so4_out, local_pow_so4_deriv_out, iconf)

module module_compute_pow_so4 

contains
 
subroutine compute_pow_so4(i_start_at, i_final_at, d_n_neigh, d_kind_neigh, iconf)

  use module_kind_variables, ONLY: double
#ifdef MLD_NDM
  use gen_com_m, ONLY: A2cm, lperiod
  use gen_com_m_ml, ONLY: imm, bg, at
  use tab_imm_m_ml, ONLY: xp
#else
  use ondm_gen_com_m, ONLY: imm, A2cm, lperiod, bg, at, indi2
  use ondm_tab_imm_m, ONLY: iwmax2, xp
#endif
  use angular_functions   ! , only : Umm
  use derived_types, only: config_real, config_desc
  use ml_in_ndm_module, ONLY: rangml, debug, imm_neigh,  rangml, j_max, jj_max, pow_so4_dim, weighted, weighted_3ch, &
                              desc_forces, linvisible, lmask 
  use module_bispectrum_so4, only: cmm, cmm2, cmm3, class_mml, &
                                   czero, fcut, dfcut, fcut_w, dfcut_w, fcut_w_3ch, dfcut_w_3ch
  use module_neigh_local, only: r_cut, r_cut_width, r_cut_in, r_cut_width_in, r_cut_pair_in, fcut_rij

#ifdef MLD_NDM
  use notperiod_mod
#else
  use ondm_transform_coord, only: ondm_notperiod
#endif
  implicit none

  integer, intent(in)  :: i_start_at, i_final_at
  integer, dimension(imm), intent(out)   :: d_n_neigh
  integer, dimension(imm, imm_neigh), intent(out)    :: d_kind_neigh
  ! double precision,dimension(0:jj_max,imm),intent(out) :: local_pow_so4_out
  ! double precision,dimension(0:jj_max,imm,0:imm_neigh,3), intent(out) :: local_pow_so4_deriv_out
  integer, optional    :: iconf

  logical  :: small
  double complex, dimension(-jj_max:jj_max, -jj_max:jj_max, 0:jj_max)  :: Umm
  double complex, dimension(-jj_max:jj_max, -jj_max:jj_max, 0:jj_max, 3)     :: dUmm
  ! cmm derivatives for the componenets for 4D spherical functions.
  ! cmm (m1,m2,j), with m1,2=-j,j, j=0,2*j_max
  ! please note that for derivatives wer have the components cmm(m1,m2,j,ia_n,3) with ia_n  running from 0 to max of neighbours of a central atom.
  ! ia_n =  0 is for central atom
  ! ia_n != 0 is for any other atom not central
  ! double complex, allocatable, dimension(:,:,:) :: cmm
  ! double complex, allocatable, dimension(:,:,:,:,:) :: dcmm

  real(double), dimension(:, :), allocatable   :: xpnp
  real(double), dimension(3) :: dxp_ji, ds
  integer  :: ia, ja, ia_n, j, iw, iw1, iw2, icnt, icnt2, icnt3

  double complex :: ctmp, ctmp2, ctmp3, sum_bi, sum_bi_w, sum_bi_w_3ch
  double complex, dimension(:), allocatable    :: dsum_bi_ia, dsum_bi_ia_w, dsum_bi_ia_w_3ch
  double complex, dimension(:), pointer  :: read_vecall

  integer  :: n_axpy
  integer, parameter   :: incx = 1, incy = 1

  integer  :: l, m1, m2
  double precision     :: r2_ji, r_ji, etmp, etmp2, etmp3
  double precision     :: factor_ia, factor_ja, factor_ia_3ch, factor_ja_3ch
  real(double) :: fcut_in, dfcut_in, fcut_out, dfcut_out 
  real(double) :: local_r_cut_in

  logical  :: desc_forces_local


  if ((i_start_at == 0) .and. (i_final_at == 0)) then
    d_n_neigh(:) = 0
    d_kind_neigh(:, :) = 0
    config_desc(iconf)%energy(:, :) = 0.d0
    ! config_desc(iconf)%force(:,:,:,:)=0.d0
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

  d_n_neigh(:) = 0
  d_kind_neigh(:, :) = 0
  config_desc(iconf)%energy(:, :) = 0.d0
  if (desc_forces_local) then
    config_desc(iconf)%force(:, :, :, :) = 0.d0
    if (allocated(dsum_bi_ia)) deallocate (dsum_bi_ia); allocate (dsum_bi_ia(imm_neigh*3))
    if (weighted) then
      if (allocated(dsum_bi_ia_w)) deallocate (dsum_bi_ia_w); allocate (dsum_bi_ia_w(imm_neigh*3))
    end if
    if (weighted_3ch) then
      if (allocated(dsum_bi_ia_w_3ch)) deallocate (dsum_bi_ia_w_3ch); allocate (dsum_bi_ia_w_3ch(imm_neigh*3))
    end if
  end if
#ifdef MLD_NDM
iw2=0
!!$    if (i_start_at == 1) iw2 = 0
!!$    if (i_start_at > 1) iw2 = iwmax2(i_start_at - 1)
#else
  if (i_start_at == 1) iw2 = 0
  if (i_start_at > 1) iw2 = iwmax2(i_start_at - 1)
#endif

  do ja = i_start_at, i_final_at
    if (lmask) then 
      if (.not.config_desc(iconf)%amask(ja)) cycle 
    end if 
    if (weighted) then
      factor_ja = config_real(iconf)%weight_per_type(config_real(iconf)%itype(ja))
      if (weighted_3ch) factor_ja_3ch = config_real(iconf)%weight_per_type_3ch(config_real(iconf)%itype(ja))
    else
      factor_ja = 1.d0
    end if

    if (linvisible .and. weighted) then
      if (config_real(iconf)%invisible_per_type(config_real(iconf)%itype(ja))) cycle
    end if
    !begin small box or not 1/
#ifdef MLD_NDM
!!$    if (small) then
#else
    if (small) then
#endif
      iw1 = 1
      iw2 = config_real(iconf)%n_neigh(ja)
#ifdef MLD_NDM
!!$    else
!!$      iw1 = iw2 + 1
!!$      iw2 = iwmax2(ja)
!!$    end if
#else
    else
      iw1 = iw2 + 1
      iw2 = iwmax2(ja)
    end if
#endif
    !end   small box or not 1/

    if (debug) then
      if (mod(ja - 1, 10) == 0) then
        if (rangml == 0) write (6, '("in bso4 i_start_at i_final_at  ja:  ",3i9)') i_start_at, i_final_at, ja
      end if
    end if

    ! INIT cmm, cmm2, cmm3
    ! cmm(:,:,:) = 0.d0
    do j = 0, jj_max; do m1 = -j, j, 2; do m2 = -j, j, 2
        cmm(m2, m1, j) = 0.d0
      end do; end do; end do
    if (weighted) then
      ! cmm2(:,:,:)=0.d0
      do j = 0, jj_max; do m1 = -j, j, 2; do m2 = -j, j, 2
          cmm2(m2, m1, j) = 0.d0
        end do; end do; end do

      if (weighted_3ch) then
        ! cmm3(:,:,:)=0.d0
        do j = 0, jj_max; do m1 = -j, j, 2; do m2 = -j, j, 2
            cmm3(m2, m1, j) = 0.d0
          end do; end do; end do
      end if
    end if

    do j = 0, jj_max
      do m1 = -j, j, 2
        cmm(m1, m1, j) = 1.d0
        if (weighted) then
          cmm2(m1, m1, j) = factor_ja                      ! 1.d0
          if (weighted_3ch) cmm3(m1, m1, j) = factor_ja_3ch                         ! 1.d0
        end if
      end do
    end do

    !DCOS    !set-up initialization of some variables for descriptors ...
    !DCOS    cmm(:,:,:) = 0.d0
    !DCOS    do j = 0,jj_max
    !DCOS      do m1=-j,j,2
    !DCOS        cmm(m1,m1,j) = 1.d0
    !DCOS      end do
    !DCOS    end do
    !DCOS    if (desc_forces_local) dcmm(:,:,:,:,:) = 0.d0

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
      !end small box or not 2/

      if (weighted) then
        if (linvisible .and. weighted) then
          if (config_real(iconf)%invisible_per_type(config_real(iconf)%itype(ia))) cycle
        end if
        factor_ia = config_real(iconf)%weight_per_type(config_real(iconf)%itype(ia))                       ! /factor_weight_mass
        if (weighted_3ch) factor_ia_3ch = config_real(iconf)%weight_per_type_3ch(config_real(iconf)%itype(ia))
      else
        factor_ia = 1.d0
      end if

#ifdef MLD_NDM
!!$      if (small) then
#else
      if (small) then
#endif
        r_ji = config_real(iconf)%r_ij(ja, iw)
        r2_ji = r_ji**2
        dxp_ji(:) = config_real(iconf)%u_ij(ja, iw, :)
#ifdef MLD_NDM
!!$      else
!!$        dxp_ji(1:3) = xpnp(1:3, ia) - xpnp(1:3, ja)
!!$        ds(:) = MatMul(dxp_ji(:), bg(:, :))
!!$        WHERE ((ds .GT. 0.5d0) .OR. (ds .LT. -0.5d0))
!!$          ds(1:3) = ds(1:3) - Dble(Nint(ds(1:3)))
!!$        END WHERE
!!$        dxp_ji(:) = MatMul(at(:, :), ds(:))
!!$        dxp_ji(:) = dxp_ji(:)/A2cm
!!$        r2_ji = Sum(dxp_ji(1:3)**2)
!!$        r_ji = dsqrt(r2_ji)
!!$      end if
#else
      else
        dxp_ji(1:3) = xpnp(1:3, ia) - xpnp(1:3, ja)
        ds(:) = MatMul(dxp_ji(:), bg(:, :))
        WHERE ((ds .GT. 0.5d0) .OR. (ds .LT. -0.5d0))
          ds(1:3) = ds(1:3) - Dble(Nint(ds(1:3)))
        END WHERE
        dxp_ji(:) = MatMul(at(:, :), ds(:))
        dxp_ji(:) = dxp_ji(:)/A2cm
        r2_ji = Sum(dxp_ji(1:3)**2)
        r_ji = dsqrt(r2_ji)
      end if
#endif

      if (r_ji >= r_cut) cycle
      ia_n = ia_n + 1

      !$! fcut = 0.5d0*(cos(one_pi*r_ji/r_cut) + 1.d0)
      !$! dfcut = -0.5d0*one_pi/r_cut*sin(one_pi*r_ji/r_cut)

      !$! begin the new fcut .....
      call fcut_rij(2, r_ji, r_cut, r_cut_width, desc_forces_local, fcut_out, dfcut_out)
      ! Pair-specific inner cutoff
      local_r_cut_in = r_cut_pair_in(config_real(iconf)%itype_db(ja), config_real(iconf)%itype_db(ia))
      if (local_r_cut_in > 0 ) then 
        call fcut_rij(3, r_ji, local_r_cut_in + r_cut_width_in, r_cut_width_in, desc_forces_local, fcut_in, dfcut_in)
      else 
        fcut_in = 0.d0 
        dfcut_in = 0.d0 
      end if  
      fcut = (1.d0 - fcut_in) * fcut_out 
      if (desc_forces_local) dfcut =-dfcut_in * fcut_out +  (1.d0 - fcut_in) * dfcut_out
      !$! end the new fcut .....



      if (weighted) then
        fcut_w = fcut*factor_ia
        dfcut_w = dfcut*factor_ia
        if (weighted_3ch) then
          fcut_w_3ch = fcut*factor_ia_3ch
          dfcut_w_3ch = dfcut*factor_ia_3ch
        end if
      end if

      call spherical_4d(ia_n, dxp_ji, r_ji, Umm, dUmm)
      ! cmm2(:,:,:)= cmm2(:,:,:)+Umm(:,:,:)*fcut_w
      do j = 0, jj_max; do m1 = -j, j, 2; do m2 = -j, j, 2
          cmm(m2, m1, j) = cmm(m2, m1, j) + Umm(m2, m1, j)*fcut
        end do; end do; end do
      if (weighted) then
        ! cmm2(:,:,:)= cmm2(:,:,:)+Umm(:,:,:)*fcut_w
        do j = 0, jj_max; do m1 = -j, j, 2; do m2 = -j, j, 2
            cmm2(m2, m1, j) = cmm2(m2, m1, j) + Umm(m2, m1, j)*fcut_w
          end do; end do; end do
        if (weighted_3ch) then
          ! cmm3(:,:,:)= cmm3(:,:,:)+Umm(:,:,:)*fcut_w_3ch
          do j = 0, jj_max; do m1 = -j, j, 2; do m2 = -j, j, 2
              cmm3(m2, m1, j) = cmm3(m2, m1, j) + Umm(m2, m1, j)*fcut_w_3ch
            end do; end do; end do
        end if
      end if
      ! cmm(:,:,:)= cmm(:,:,:)+Umm(:,:,:)*fcut
      ! if (desc_forces_local) then
      !  do ix=1,3
      !    dcmm(:,:,:,ia_n,ix) =  (dfcut*dxp_ji(ix)*Umm(:,:,:)/r_ji + fcut*dUmm(:,:,:,ix))
      !  end do
      ! end if
      d_kind_neigh(ja, ia_n) = ia
    end do                  ! iw
    d_n_neigh(ja) = ia_n

    !DCOS if (desc_forces_local) then
    !DCOS   dcmm(:,:,:,0,:) = - sum(dcmm(:,:,:,:,:), dim=4)
    !DCOS end if

    icnt = 0
    n_axpy = 3*ia_n
    do l = 0, jj_max
      icnt = icnt + 1
      etmp = 0.d0
      etmp2 = 0.d0
      etmp3 = 0.d0
      if (desc_forces_local) then
        dsum_bi_ia(1:3*ia_n) = czero
        if (weighted) then
          dsum_bi_ia_w(1:3*ia_n) = czero
          if (weighted_3ch) then
            dsum_bi_ia_w_3ch(1:3*ia_n) = czero
          end if
        end if
      end if
      do m1 = -l, l, 2
        do m2 = -l, l, 2
          ctmp = cmm(m2, m1, l)
          ! config_desc(iconf)%energy(icnt,ja) = config_desc(iconf)%energy(icnt,ja) + real(conjg(ctmp) * ctmp)
          etmp = etmp + real(conjg(ctmp)*ctmp)
          if (weighted) then
            ctmp2 = cmm2(m2, m1, l)
            etmp2 = etmp2 + real(conjg(ctmp2)*ctmp2)
            if (weighted_3ch) then
              ctmp3 = cmm3(m2, m1, l)
              etmp3 = etmp3 + real(conjg(ctmp3)*ctmp3)
            end if
          end if
          if (desc_forces_local) then
            sum_bi = 2.d0*ctmp
            read_vecall => class_mml(m2, m1, l)%vecall
            ! do ia=1,ia_n
            !  config_desc(iconf)%force(icnt,ja,ia,1:3) = config_desc(iconf)%force(icnt,ja,ia,1:3) + conjg(dcmm(m2,m1,l,ia,1:3))*sum_bi
            ! end do !ia
            call zaxpy(n_axpy, sum_bi, read_vecall, incx, dsum_bi_ia, incy)
            nullify (read_vecall)
            if (weighted) then
              sum_bi_w = 2.d0*ctmp2
              read_vecall => class_mml(m2, m1, l)%vecall_w
              call zaxpy(n_axpy, sum_bi_w, read_vecall, incx, dsum_bi_ia_w, incy)
              nullify (read_vecall)
            end if
            if (weighted_3ch) then
              sum_bi_w_3ch = 2.d0*ctmp3
              read_vecall => class_mml(m2, m1, l)%vecall_w_3ch
              call zaxpy(n_axpy, sum_bi_w_3ch, read_vecall, incx, dsum_bi_ia_w_3ch, incy)
              nullify (read_vecall)
            end if
          end if                  ! desc_forces_local
        end do                  ! m2
      end do                  ! m1

      ! fix energy descriptor
      config_desc(iconf)%energy(icnt, ja) = etmp
      if (weighted) then
        icnt2 = pow_so4_dim + icnt
        icnt3 = pow_so4_dim + icnt2
        config_desc(iconf)%energy(icnt2, ja) = etmp2
        if (weighted_3ch) then
          config_desc(iconf)%energy(icnt3, ja) = etmp3
        end if
      end if
      ! fix forces descriptor
      if (desc_forces_local) then
        do ia = 1, ia_n
          config_desc(iconf)%force(icnt, ja, ia, 1) = real(dsum_bi_ia(3*ia - 2), kind=kind(1.d0))
          config_desc(iconf)%force(icnt, ja, ia, 2) = real(dsum_bi_ia(3*ia - 1), kind=kind(1.d0))
          config_desc(iconf)%force(icnt, ja, ia, 3) = real(dsum_bi_ia(3*ia), kind=kind(1.d0))
        end do
        if (weighted) then
          do ia = 1, ia_n
            config_desc(iconf)%force(icnt2, ja, ia, 1) = real(dsum_bi_ia_w(3*ia - 2), kind=kind(1.d0))
            config_desc(iconf)%force(icnt2, ja, ia, 2) = real(dsum_bi_ia_w(3*ia - 1), kind=kind(1.d0))
            config_desc(iconf)%force(icnt2, ja, ia, 3) = real(dsum_bi_ia_w(3*ia), kind=kind(1.d0))
          end do
          if (weighted_3ch) then
            do ia = 1, ia_n
              config_desc(iconf)%force(icnt3, ja, ia, 1) = real(dsum_bi_ia_w_3ch(3*ia - 2), kind=kind(1.d0))
              config_desc(iconf)%force(icnt3, ja, ia, 2) = real(dsum_bi_ia_w_3ch(3*ia - 1), kind=kind(1.d0))
              config_desc(iconf)%force(icnt3, ja, ia, 3) = real(dsum_bi_ia_w_3ch(3*ia), kind=kind(1.d0))
            end do
          end if
        end if
      end if

    end do                  ! l

    ! if (desc_forces_local)  config_desc(iconf)%force(:,ja,0,:) = -sum(config_desc(iconf)%force(:,ja,1:ia_n,:), dim=2)
    ! if (desc_forces_local)  then
    !   write (*,*) shape(config_desc(iconf)%force(:,ja,:,:))
    !   config_desc(iconf)%force(1:pow_so4_dim,ja,0,1:3) = -sum(config_desc(iconf)%force(1:pow_so4_dim,ja,1:ia_n,1:3), dim=2)
    ! end if
    if (desc_forces_local) then
    do ia = 1, ia_n
      config_desc(iconf)%force(:, ja, 0, :) = config_desc(iconf)%force(:, ja, 0, :) - config_desc(iconf)%force(:, ja, ia, :)
    end do
    end if

    ! config_desc(iconf)%energy(1:pow_so4_dim,ja) = config_desc(iconf)%energy(1:pow_so4_dim,ja)
    if (weighted) then
      config_desc(iconf)%energy(pow_so4_dim + 1:2*pow_so4_dim, ja) = config_desc(iconf)%energy(pow_so4_dim + 1:2*pow_so4_dim, ja) ! *factor_ja
      if (weighted_3ch) config_desc(iconf)%energy(2*pow_so4_dim + 1:3*pow_so4_dim, ja) = config_desc(iconf)%energy(2*pow_so4_dim + 1:3*pow_so4_dim, ja)    ! *factor_ja_3ch
      if (desc_forces_local) then
        config_desc(iconf)%force(pow_so4_dim + 1:2*pow_so4_dim, ja, :, :) = config_desc(iconf)%force(pow_so4_dim + 1:2*pow_so4_dim, ja, :, :)                ! *factor_ja
        if (weighted_3ch) config_desc(iconf)%force(2*pow_so4_dim + 1:3*pow_so4_dim, ja, :, :) = config_desc(iconf)%force(2*pow_so4_dim + 1:3*pow_so4_dim, ja, :, :)                   ! *factor_ja_3ch
      end if
    end if

  end do                  ! ja

  deallocate (xpnp)

end subroutine compute_pow_so4
end module module_compute_pow_so4 




subroutine init_pow_so4()
  use ml_in_ndm_module, only: jj_max, imm_neigh, weighted, weighted_3ch
  use module_bispectrum_so4, only: class_mml, cmm, cmm2, cmm3
  implicit none

  integer  :: j, m1, m2


  if (allocated(class_mml)) deallocate (class_mml)
  allocate (class_mml(-jj_max:jj_max, -jj_max:jj_max, 0:jj_max))
  do j = 0, jj_max
    do m1 = -j, j, 2
      do m2 = -j, j, 2
        if (allocated(class_mml(m2, m1, j)%vecall)) deallocate (class_mml(m2, m1, j)%vecall)
        allocate (class_mml(m2, m1, j)%vecall(1:3*imm_neigh))
      end do
    end do
  end do

  if (allocated(cmm)) deallocate (cmm)
  allocate (cmm(-jj_max:jj_max, -jj_max:jj_max, 0:jj_max))
  if (weighted) then
    if (allocated(cmm2)) deallocate (cmm2)
    allocate (cmm2(-jj_max:jj_max, -jj_max:jj_max, 0:jj_max))
    do j = 0, jj_max
      do m1 = -j, j, 2
        do m2 = -j, j, 2
          if (allocated(class_mml(m2, m1, j)%vecall_w)) deallocate (class_mml(m2, m1, j)%vecall_w)
          allocate (class_mml(m2, m1, j)%vecall_w(1:3*imm_neigh))
        end do
      end do
    end do

    if (weighted_3ch) then
      if (allocated(cmm3)) deallocate (cmm3)
      allocate (cmm3(-jj_max:jj_max, -jj_max:jj_max, 0:jj_max))
      do j = 0, jj_max
        do m1 = -j, j, 2
          do m2 = -j, j, 2
            if (allocated(class_mml(m2, m1, j)%vecall_w_3ch)) deallocate (class_mml(m2, m1, j)%vecall_w_3ch)
            allocate (class_mml(m2, m1, j)%vecall_w_3ch(1:3*imm_neigh))
          end do
        end do
      end do
    end if                  ! weighted_3ch
  end if                  ! weighted

end subroutine init_pow_so4
