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

module angular_functions

  use ml_in_ndm_module, only: j_max, jj_max, rangml, r_factorial

  implicit none

  public   :: spherical_harm, grad_spherical_harm
contains

  subroutine spherical_4d(ia_n, x, rr, Umm, dUmm)
    ! Input:
    !    x(3) -> x, y, z ;
    !    rr=sqrt(x**2 + y**2 + z**2)
    ! Output:
    !    Spherical 4D and derivatives
    !    Umm, dUmm
    use module_kind_variables, only: kind_double
    use ml_in_ndm_module, only: jj_max, weighted, weighted_3ch
    use module_bispectrum_so4, only: inv_r0, r0, class_mml, &
                                     fcut, dfcut, &
                                     fcut_w, dfcut_w, &
                                     fcut_w_3ch, dfcut_w_3ch
    use module_neigh_local, only: r_cut 


    integer, intent(in)  :: ia_n
    real(kind_double), dimension(3), intent(in)  :: x
    real(kind_double), intent(in)    :: rr
    double complex, dimension(-jj_max:jj_max, -jj_max:jj_max, 0:jj_max), intent(out) :: Umm
    double complex, dimension(-jj_max:jj_max, -jj_max:jj_max, 0:jj_max, 1:3), intent(out)  :: dUmm
    !double complex,dimension(-jj_max:jj_max,-jj_max:jj_max,0:jj_max,3) :: dUmm

    integer  :: j, m1, m2
    real(kind_double)    :: z0, th0, l0, l0i, dz, dil0, tmp_tan_i, tmp_sin_i, &
                            tmp_xx, tmp_yy, tmp_zz, &
                            tmp_xx_w, tmp_yy_w, tmp_zz_w, &
                            tmp_xx_w_3ch, tmp_yy_w_3ch, tmp_zz_w_3ch
    real(kind_double), dimension(3)  :: drr
    double complex :: cplx_l0i, cplx_i, z_plus, z_minus, x_plus, x_minus, tjm1, tjm2, &
      tmp_m2p1_m1m1_jm1, tmp_m2m1_m1m1_jm1, tmp_m2m1j, cttt, &
      tmp_m2m1_m1p1_jm1, tmp_m2p1_m1p1_jm1


    double complex, dimension(3)     :: dz_plus, dz_minus, dx_plus, dx_minus, &
      dtmp_m2p1_m1m1_jm1, dtmp_m2m1_m1m1_jm1, dtmp_m2m1j

    cplx_i = cmplx(0.d0, 1.d0, kind=kind(1.d0))
    r0 = r_cut/inv_r0
    th0 = rr/r0
    tmp_sin_i = 1.d0/dsin(th0)
    tmp_tan_i = 1.d0/dtan(th0)
    l0 = rr*tmp_sin_i
    l0i = 1.d0/l0
    z0 = rr*tmp_tan_i

    ! z_+/- is cos(theta0) +/- i sin(theta0) * cos(theta)
    z_plus = cmplx(z0, x(3), kind=kind(1.d0))*l0i
    z_minus = cmplx(z0, -x(3), kind=kind(1.d0))*l0i
    ! x +/- is sin(theta) * sin(theta0) * exp(+/- i \phi)
    x_plus = cmplx(x(1), x(2), kind=kind(1.d0))*l0i
    x_minus = cmplx(x(1), -x(2), kind=kind(1.d0))*l0i

    !Here are the derivatives ...
    drr(1:3) = x(1:3)/rr

    !error dz  = (1.d0/tan(th0) - th0/tan(th0)**2)
    dz = (tmp_tan_i - th0*tmp_sin_i**2)
    dil0 = (dcos(th0)/r0 - l0i)/rr
    cplx_l0i = cmplx(0.d0, l0i, kind=kind(1.d0))

    dz_plus(1:3) = (cmplx(z0, x(3), kind=kind(1.d0))*dil0 + dz*l0i)*drr(1:3)
    dz_plus(3) = dz_plus(3) + cplx_l0i

    dz_minus(1:3) = (cmplx(z0, -x(3), kind=kind(1.d0))*dil0 + dz*l0i)*drr(1:3)
    dz_minus(3) = dz_minus(3) - cplx_l0i


    dx_plus(1:3) = cmplx(x(1), x(2), kind=kind(1.d0))*dil0*drr(1:3)
    dx_plus(1) = dx_plus(1) + cmplx(l0i, 0.d0, kind=kind(1.d0))
    dx_plus(2) = dx_plus(2) + cplx_l0i

    dx_minus(1:3) = cmplx(x(1), -x(2), kind=kind(1.d0))*dil0*drr(1:3)
    dx_minus(1) = dx_minus(1) + cmplx(l0i, 0.d0, kind=kind(1.d0))
    dx_minus(2) = dx_minus(2) - cplx_l0i

    do j = 1, jj_max; do m1 = -j, j, 2; do m2 = -j, j, 2
        Umm(m2, m1, j) = (0.d0, 0.d0)
      end do; end do; end do
    Umm(0, 0, 0) = (1.d0, 0.d0)
    !TODObungabunga
    dUmm(0, 0, 0, :) = (0.d0, 0.d0)
    do j = 1, jj_max
      do m1 = -j, j, 2
        Umm(m1, m1, j) = (1.d0, 0.d0)
      end do
    end do
    tmp_xx = dfcut*drr(1)
    tmp_yy = dfcut*drr(2)
    tmp_zz = dfcut*drr(3)
    class_mml(0, 0, 0)%vecall(3*ia_n - 2) = tmp_xx
    class_mml(0, 0, 0)%vecall(3*ia_n - 1) = tmp_yy
    class_mml(0, 0, 0)%vecall(3*ia_n) = tmp_zz

    if (weighted) then
      tmp_xx_w = dfcut_w*drr(1)
      tmp_yy_w = dfcut_w*drr(2)
      tmp_zz_w = dfcut_w*drr(3)
      class_mml(0, 0, 0)%vecall_w(3*ia_n - 2) = tmp_xx_w
      class_mml(0, 0, 0)%vecall_w(3*ia_n - 1) = tmp_yy_w
      class_mml(0, 0, 0)%vecall_w(3*ia_n) = tmp_zz_w
      if (weighted_3ch) then
        tmp_xx_w_3ch = dfcut_w_3ch*drr(1)
        tmp_yy_w_3ch = dfcut_w_3ch*drr(2)
        tmp_zz_w_3ch = dfcut_w_3ch*drr(3)
        class_mml(0, 0, 0)%vecall_w_3ch(3*ia_n - 2) = tmp_xx_w_3ch
        class_mml(0, 0, 0)%vecall_w_3ch(3*ia_n - 1) = tmp_yy_w_3ch
        class_mml(0, 0, 0)%vecall_w_3ch(3*ia_n) = tmp_zz_w_3ch
      end if
    end if

    do j = 1, jj_max
      do m1 = -j, j, 2
        do m2 = -j, j, 2
          if (m1 .eq. j) then
            if (m2 .eq. -j) then
              tjm2 = cplx_i*dsqrt(dble(j - m2)/dble(j + m1))
              tmp_m2p1_m1m1_jm1 = Umm(m2 + 1, m1 - 1, j - 1)
              dtmp_m2p1_m1m1_jm1(:) = dUmm(m2 + 1, m1 - 1, j - 1, :)
              ! Umm(m2,m1,j) = - tjm2*x_plus*Umm(m2+1,m1-1,j-1)
              ! dUmm(m2,m1,j,:) = - tjm2*( dx_plus(:)*Umm(m2+1,m1-1,j-1)+x_plus*dUmm(m2+1,m1-1,j-1,:) )
              tmp_m2m1j = -tjm2*x_plus*tmp_m2p1_m1m1_jm1
              dtmp_m2m1j(:) = -tjm2*(dx_plus(:)*tmp_m2p1_m1m1_jm1 + x_plus*dtmp_m2p1_m1m1_jm1(:))
            else if (m2 .eq. j) then
              tjm1 = dsqrt(dble(j + m2)/dble(j + m1))
              tmp_m2m1_m1m1_jm1 = Umm(m2 - 1, m1 - 1, j - 1)
              dtmp_m2m1_m1m1_jm1(:) = dUmm(m2 - 1, m1 - 1, j - 1, :)
              ! Umm(m2,m1,j) = tjm1*z_minus*Umm(m2-1,m1-1,j-1)
              ! dUmm(m2,m1,j,:) = tjm1*(dz_minus(:)*Umm(m2-1,m1-1,j-1)+z_minus*dUmm(m2-1,m1-1,j-1,:))
              tmp_m2m1j = tjm1*z_minus*tmp_m2m1_m1m1_jm1
              dtmp_m2m1j(:) = tjm1*(dz_minus(:)*tmp_m2m1_m1m1_jm1 + z_minus*dtmp_m2m1_m1m1_jm1(:))
            else
              tjm1 = dsqrt(dble(j + m2)/dble(j + m1))
              tjm2 = cplx_i*dsqrt(dble(j - m2)/dble(j + m1))
              tmp_m2m1_m1m1_jm1 = Umm(m2 - 1, m1 - 1, j - 1)
              dtmp_m2m1_m1m1_jm1(:) = dUmm(m2 - 1, m1 - 1, j - 1, :)
              tmp_m2p1_m1m1_jm1 = Umm(m2 + 1, m1 - 1, j - 1)
              dtmp_m2p1_m1m1_jm1(:) = dUmm(m2 + 1, m1 - 1, j - 1, :)
              ! Umm(m2,m1,j) = tjm1*z_minus*Umm(m2-1,m1-1,j-1) - tjm2*x_plus*Umm(m2+1,m1-1,j-1)
              ! dUmm(m2,m1,j,:) = tjm1*(dz_minus(:)*Umm(m2-1,m1-1,j-1) + z_minus*dUmm(m2-1,m1-1,j-1,:))   &
              !                  -tjm2*(dx_plus(:)*Umm(m2+1,m1-1,j-1)+x_plus*dUmm(m2+1,m1-1,j-1,:))
              tmp_m2m1j = tjm1*z_minus*tmp_m2m1_m1m1_jm1 - tjm2*x_plus*tmp_m2p1_m1m1_jm1
              dtmp_m2m1j(:) = tjm1*(dz_minus(:)*tmp_m2m1_m1m1_jm1 + z_minus*dtmp_m2m1_m1m1_jm1(:)) &
                              - tjm2*(dx_plus(:)*tmp_m2p1_m1m1_jm1 + x_plus*dtmp_m2p1_m1m1_jm1(:))
            end if
          else if (m1 .eq. -j) then
            if (m2 .eq. j) then
              tjm1 = cplx_i*dsqrt(dble(j + m2)/dble(j - m1))
              tmp_m2m1_m1p1_jm1 = Umm(m2 - 1, m1 + 1, j - 1)
              ! Umm(m2,m1,j) =  - tjm1*x_minus*Umm(m2-1,m1+1,j-1)
              ! dUmm(m2,m1,j,:) =  - tjm1*(dx_minus(:)*Umm(m2-1,m1+1,j-1)+x_minus*dUmm(m2-1,m1+1,j-1,:))
              tmp_m2m1j = -tjm1*x_minus*tmp_m2m1_m1p1_jm1
              dtmp_m2m1j(:) = -tjm1*(dx_minus(:)*tmp_m2m1_m1p1_jm1 + x_minus*dUmm(m2 - 1, m1 + 1, j - 1, :))
            else if (m2 .eq. -j) then
              tjm2 = dsqrt(dble(j - m2)/dble(j - m1))
              tmp_m2p1_m1p1_jm1 = Umm(m2 + 1, m1 + 1, j - 1)
              ! Umm(m2,m1,j) = tjm2*z_plus*Umm(m2+1,m1+1,j-1)
              ! dUmm(m2,m1,j,:) = tjm2*(dz_plus(:)*Umm(m2+1,m1+1,j-1)+z_plus*dUmm(m2+1,m1+1,j-1,:))
              tmp_m2m1j = tjm2*z_plus*tmp_m2p1_m1p1_jm1
              ! write (*,*) 'tjm2', tjm2
              ! write (*,*) 'tmp_m2p1_m1p1_jm1', tmp_m2p1_m1p1_jm1
              ! write (*,*) 'dz_plus', dz_plus(:)
              ! write (*,*) 'z_plus', z_plus,m2,m1,j
              ! write (*,*) 'dUmm', dUmm(m2+1,m1+1,j-1,:)
              dtmp_m2m1j(:) = tjm2*(dz_plus(:)*tmp_m2p1_m1p1_jm1 + z_plus*dUmm(m2 + 1, m1 + 1, j - 1, :))
            else

              tjm1 = cplx_i*dsqrt(dble(j + m2)/dble(j - m1))
              tjm2 = dsqrt(dble(j - m2)/dble(j - m1))
              tmp_m2m1_m1p1_jm1 = Umm(m2 - 1, m1 + 1, j - 1)
              tmp_m2p1_m1p1_jm1 = Umm(m2 + 1, m1 + 1, j - 1)
              ! Umm(m2,m1,j) = tjm2*z_plus*Umm(m2+1,m1+1,j-1) - tjm1*x_minus*Umm(m2-1,m1+1,j-1)
              ! dUmm(m2,m1,j,:) = tjm2*(dz_plus(:)*Umm(m2+1,m1+1,j-1)+z_plus*dUmm(m2+1,m1+1,j-1,:)) &
              !                 - tjm1*(dx_minus(:)*Umm(m2-1,m1+1,j-1)+x_minus*dUmm(m2-1,m1+1,j-1,:))
              tmp_m2m1j = tjm2*z_plus*tmp_m2p1_m1p1_jm1 - tjm1*x_minus*tmp_m2m1_m1p1_jm1
              dtmp_m2m1j(:) = tjm2*(dz_plus(:)*tmp_m2p1_m1p1_jm1 + z_plus*dUmm(m2 + 1, m1 + 1, j - 1, :)) &
                              - tjm1*(dx_minus(:)*tmp_m2m1_m1p1_jm1 + x_minus*dUmm(m2 - 1, m1 + 1, j - 1, :))

            end if
          else
            if (m2 .eq. -j) then
              tjm2 = cplx_i*dsqrt(dble(j - m2)/dble(j + m1))
              tmp_m2p1_m1m1_jm1 = Umm(m2 + 1, m1 - 1, j - 1)
              dtmp_m2p1_m1m1_jm1(:) = dUmm(m2 + 1, m1 - 1, j - 1, :)
              ! Umm(m2,m1,j) = - tjm2*x_plus*Umm(m2+1,m1-1,j-1)
              ! dUmm(m2,m1,j,:) = - tjm2*( dx_plus(:)*Umm(m2+1,m1-1,j-1)+x_plus*dUmm(m2+1,m1-1,j-1,:) )
              tmp_m2m1j = -tjm2*x_plus*tmp_m2p1_m1m1_jm1
              dtmp_m2m1j(:) = -tjm2*(dx_plus(:)*tmp_m2p1_m1m1_jm1 + x_plus*dtmp_m2p1_m1m1_jm1(:))
            else if (m2 .eq. j) then
              tjm1 = dsqrt(dble(j + m2)/dble(j + m1))
              tmp_m2m1_m1m1_jm1 = Umm(m2 - 1, m1 - 1, j - 1)
              dtmp_m2m1_m1m1_jm1(:) = dUmm(m2 - 1, m1 - 1, j - 1, :)
              ! Umm(m2,m1,j) = tjm1*z_minus*Umm(m2-1,m1-1,j-1)
              ! dUmm(m2,m1,j,:) = tjm1*(dz_minus(:)*Umm(m2-1,m1-1,j-1)+z_minus*dUmm(m2-1,m1-1,j-1,:))
              tmp_m2m1j = tjm1*z_minus*tmp_m2m1_m1m1_jm1
              dtmp_m2m1j(:) = tjm1*(dz_minus(:)*tmp_m2m1_m1m1_jm1 + z_minus*dtmp_m2m1_m1m1_jm1(:))
            else
              tjm1 = dsqrt(dble(j + m2)/dble(j + m1))
              tjm2 = cplx_i*dsqrt(dble(j - m2)/dble(j + m1))
              tmp_m2m1_m1m1_jm1 = Umm(m2 - 1, m1 - 1, j - 1)
              dtmp_m2m1_m1m1_jm1(:) = dUmm(m2 - 1, m1 - 1, j - 1, :)
              tmp_m2p1_m1m1_jm1 = Umm(m2 + 1, m1 - 1, j - 1)
              dtmp_m2p1_m1m1_jm1(:) = dUmm(m2 + 1, m1 - 1, j - 1, :)
              ! Umm(m2,m1,j) = tjm1*z_minus*Umm(m2-1,m1-1,j-1) - tjm2*x_plus*Umm(m2+1,m1-1,j-1)
              ! dUmm(m2,m1,j,:) = tjm1*(dz_minus(:)*Umm(m2-1,m1-1,j-1) + z_minus*dUmm(m2-1,m1-1,j-1,:))   &
              !                  -tjm2*(dx_plus(:)*Umm(m2+1,m1-1,j-1)+x_plus*dUmm(m2+1,m1-1,j-1,:))
              tmp_m2m1j = tjm1*z_minus*tmp_m2m1_m1m1_jm1 - tjm2*x_plus*tmp_m2p1_m1m1_jm1
              dtmp_m2m1j(:) = tjm1*(dz_minus(:)*tmp_m2m1_m1m1_jm1 + z_minus*dtmp_m2m1_m1m1_jm1(:)) &
                              - tjm2*(dx_plus(:)*tmp_m2p1_m1m1_jm1 + x_plus*dtmp_m2p1_m1m1_jm1(:))
            end if
          end if

          Umm(m2, m1, j) = tmp_m2m1j
          dUmm(m2, m1, j, 1:3) = dtmp_m2m1j(1:3)
          cttt = conjg(tmp_m2m1j)
          class_mml(m2, m1, j)%vecall(3*ia_n - 2) = cttt*tmp_xx + fcut*conjg(dtmp_m2m1j(1))
          class_mml(m2, m1, j)%vecall(3*ia_n - 1) = cttt*tmp_yy + fcut*conjg(dtmp_m2m1j(2))
          class_mml(m2, m1, j)%vecall(3*ia_n) = cttt*tmp_zz + fcut*conjg(dtmp_m2m1j(3))
          if (weighted) then
            class_mml(m2, m1, j)%vecall_w(3*ia_n - 2) = cttt*tmp_xx_w + fcut_w*conjg(dtmp_m2m1j(1))
            class_mml(m2, m1, j)%vecall_w(3*ia_n - 1) = cttt*tmp_yy_w + fcut_w*conjg(dtmp_m2m1j(2))
            class_mml(m2, m1, j)%vecall_w(3*ia_n) = cttt*tmp_zz_w + fcut_w*conjg(dtmp_m2m1j(3))
            if (weighted_3ch) then
              class_mml(m2, m1, j)%vecall_w_3ch(3*ia_n - 2) = cttt*tmp_xx_w_3ch + fcut_w_3ch*conjg(dtmp_m2m1j(1))
              class_mml(m2, m1, j)%vecall_w_3ch(3*ia_n - 1) = cttt*tmp_yy_w_3ch + fcut_w_3ch*conjg(dtmp_m2m1j(2))
              class_mml(m2, m1, j)%vecall_w_3ch(3*ia_n) = cttt*tmp_zz_w_3ch + fcut_w_3ch*conjg(dtmp_m2m1j(3))
            end if
          end if

        end do
      end do
    end do
  end subroutine spherical_4d



  function solid_harm(l, m, x)
    ! From Math Wolfram World
    ! http://functions.wolfram.com/PDF/SphericalHarmonicY.pdf
    ! page 13 at the bottom
    ! Varshalovich page 133 eq.16

    integer, intent(in)  :: l, m
    double precision, dimension(3), intent(in)   :: x
    integer  :: i, j, k
    double complex :: solid_harm, ct1, ct2, ct3, tmp_ij, tmp_kk
    logical :: xy_zero, z_zero

    solid_harm = (0.d0, 0.d0)
    ct1 = dcmplx(-0.5d0*x(1), -0.5d0*x(2))
    ct2 = dcmplx(0.5d0*x(1), -0.5d0*x(2))
    ct3 = x(3)
    ! Intel __powc16i4 raises FPE on (0+0i)**n for any n.
    ! Since |x|>0, xy_zero and z_zero are mutually exclusive.
    xy_zero = (x(1) == 0.d0) .and. (x(2) == 0.d0)
    z_zero  = (x(3) == 0.d0)
    do i = 0, l
      j = i - m               ! \delta_{i-j,m}
      k = l - i - j           ! \delta_{i+j+k,l}
      !
      if ((j .ge. 0) .and. (k .ge. 0)) then
        ! ct1**i * ct2**j factor: zero unless both exponents are 0 when base is 0
        if (xy_zero) then
          if (i == 0 .and. j == 0) then; tmp_ij = (1.d0, 0.d0)
          else; cycle; end if  ! whole term is zero, skip
        else
          tmp_ij = ct1**i * ct2**j
        end if
        ! ct3**k factor
        if (z_zero) then
          if (k == 0) then; tmp_kk = (1.d0, 0.d0); else; cycle; end if
        else
          tmp_kk = ct3**k
        end if
        solid_harm = solid_harm + tmp_ij * tmp_kk &
          / (r_factorial(i)*r_factorial(j)*r_factorial(k))
      end if
    end do
    solid_harm = solid_harm*dsqrt(r_factorial(l + m)*r_factorial(l - m))
  end function solid_harm


  function spherical_harm(l, m, x)
    ! From Math Wolfram World
    ! http://functions.wolfram.com/PDF/SphericalHarmonicY.pdf
    ! page 13 at the botom
    use ml_in_ndm_module, only: one_pi

    integer, intent(in)  :: l, m
    double precision, dimension(3), intent(in)   :: x
    double complex :: spherical_harm
    double precision :: r2

    r2 = SUM(x(1:3)**2)
    spherical_harm = solid_harm(l, m, x)*dsqrt((2.d0*dble(l) + 1.d0)/(4.d0*one_pi))*r2**(-0.5d0*l)
  end function spherical_harm



  function grad_spherical_harm(l, m, x)
    ! From Math Wolfram World
    ! http://functions.wolfram.com/PDF/SphericalHarmonicY.pdf
    ! page 13 at the botom
    use ml_in_ndm_module, only: one_pi
    implicit none

    integer, intent(in)  :: l, m
    double precision, dimension(3), intent(in)   :: x
    double complex, dimension(3)     :: grad_spherical_harm
    integer  :: i, j, k
    double complex :: tmp, tmp_i, tmp_j, tmp_k, tmp_ii, tmp_jj, tmp_kk, &
      tmp_ii1, tmp_jj1, tmp_kk1
    double precision :: r2
    logical :: xy_zero, z_zero

    r2 = sum(x(1:3)**2)

    tmp_i = dcmplx(-0.5d0*x(1), -0.5d0*x(2))
    tmp_j = dcmplx(0.5d0*x(1), -0.5d0*x(2))
    tmp_k = x(3)
    ! Intel __powc16i4 raises FPE on (0+0i)**n for any n.
    ! Since |x|>0 the two cases are mutually exclusive:
    !   xy_zero: x(1)=x(2)=0 → tmp_i=tmp_j=0, but tmp_k/=0
    !   z_zero:  x(3)=0      → tmp_k=0,       but tmp_i,tmp_j/=0
    xy_zero = (x(1) == 0.d0) .and. (x(2) == 0.d0)
    z_zero  = (x(3) == 0.d0)
    grad_spherical_harm(:) = (0.d0, 0.d0)
    do i = 0, l
      j = i - m
      k = l - i - j
      ! Only guard the base that can be zero; the other is always safe.
      if (xy_zero) then
        if (i == 0) then; tmp_ii = (1.d0, 0.d0); else; tmp_ii = (0.d0, 0.d0); end if
        if (j >= 0) then
          if (j == 0) then; tmp_jj = (1.d0, 0.d0); else; tmp_jj = (0.d0, 0.d0); end if
        end if
        if (k >= 0) tmp_kk = tmp_k**k
      else if (z_zero) then
        tmp_ii = tmp_i**i
        if (j >= 0) tmp_jj = tmp_j**j
        if (k >= 0) then
          if (k == 0) then; tmp_kk = (1.d0, 0.d0); else; tmp_kk = (0.d0, 0.d0); end if
        end if
      else
        tmp_ii = tmp_i**i
        if (j >= 0) tmp_jj = tmp_j**j
        if (k >= 0) tmp_kk = tmp_k**k
      end if

      if ((i .ge. 1) .and. (j .ge. 0) .and. (k .ge. 0)) then
        if (xy_zero) then
          if (i == 1) then; tmp_ii1 = (1.d0, 0.d0); else; tmp_ii1 = (0.d0, 0.d0); end if
        else
          tmp_ii1 = tmp_i**(i - 1)
        end if
        tmp_ii1 = tmp_ii1*tmp_jj*tmp_kk/(r_factorial(i - 1)*r_factorial(j)*r_factorial(k))
        grad_spherical_harm(1) = grad_spherical_harm(1) - 0.5d0*tmp_ii1
        grad_spherical_harm(2) = grad_spherical_harm(2) - dcmplx(0.d0, 0.5d0)*tmp_ii1
      end if

      if ((i .ge. 0) .and. (j .ge. 1) .and. (k .ge. 0)) then
        if (xy_zero) then
          if (j == 1) then; tmp_jj1 = (1.d0, 0.d0); else; tmp_jj1 = (0.d0, 0.d0); end if
        else
          tmp_jj1 = tmp_j**(j - 1)
        end if
        tmp_jj1 = tmp_ii*tmp_jj1*tmp_kk/(r_factorial(i)*r_factorial(j - 1)*r_factorial(k))
        grad_spherical_harm(1) = grad_spherical_harm(1) + 0.5d0*tmp_jj1
        grad_spherical_harm(2) = grad_spherical_harm(2) - dcmplx(0.d0, 0.5d0)*tmp_jj1
      end if

      if ((i .ge. 0) .and. (j .ge. 0) .and. (k .ge. 1)) then
        if (z_zero) then
          if (k == 1) then; tmp_kk1 = (1.d0, 0.d0); else; tmp_kk1 = (0.d0, 0.d0); end if
        else
          tmp_kk1 = tmp_k**(k - 1)
        end if
        tmp = 1.d0/(r_factorial(i)*r_factorial(j)*r_factorial(k - 1))
        grad_spherical_harm(3) = grad_spherical_harm(3) + tmp_ii*tmp_jj*tmp_kk1*tmp
      end if

    end do
    grad_spherical_harm(:) = grad_spherical_harm(:)*dsqrt((2.d0*l + 1.d0)*r_factorial(l + m)*r_factorial(l - m)/(4.d0*one_pi))*r2**(-0.5d0*l)
    grad_spherical_harm(1:3) = grad_spherical_harm(1:3) - dble(l)*spherical_harm(l, m, x)*x(1:3)/r2
  end function grad_spherical_harm


end module angular_functions



function cg(j1, m1, j2, m2, j, m, ntype)
  implicit none

  integer, intent(in)  :: j1, j2, j, m1, m2, m, ntype
  double precision     :: cg
  double precision     :: wigner3j

  cg = (-1.d0)**((j1 - j2 + m)/ntype)*dsqrt(2.d0*dble(j)/dble(ntype) + 1.d0)*wigner3j(j1, m1, j2, m2, j, -m, ntype)
end function cg


function wigner3j(j1, m1, j2, m2, j, m, ntype)
  implicit none

  integer, intent(in)  :: j1, j2, j, m1, m2, m, ntype
  double precision     :: wigner3j
  integer  :: t, tmin, tmax
  double precision     :: factorial, triangle, coef, sum_coef

  coef = (-1.d0)**((j1 - j2 - m)/ntype)*dsqrt( &
         factorial((j1 + m1)/ntype)*factorial((j1 - m1)/ntype)* &
         factorial((j2 + m2)/ntype)*factorial((j2 - m2)/ntype)* &
         factorial((j + m)/ntype)*factorial((j - m)/ntype))

  triangle = dsqrt(factorial((j1 + j2 - j)/ntype)*factorial((j1 - j2 + j)/ntype)* &
                   factorial((-j1 + j2 + j)/ntype)/factorial((j1 + j2 + j)/ntype + 1))

  sum_coef = 0.d0
  tmin = max((j2 - j - m1)/ntype, (j1 + m2 - j)/ntype, 0)
  tmax = min((j1 + j2 - j)/ntype, (j1 - m1)/ntype, (j2 + m2)/ntype)
  do t = tmin, tmax
    sum_coef = sum_coef + (-1.d0)**t/( &
               factorial(t)* &
               factorial((j1 + j2 - j)/ntype - t)* &
               factorial((j1 - m1)/ntype - t)* &
               factorial((j2 + m2)/ntype - t)* &
               factorial((j - j2 + m1)/ntype + t)* &
               factorial((j - j1 - m2)/ntype + t))
  end do

  wigner3j = coef*triangle*sum_coef
end function wigner3j


function factorial(n) result(res)

  use ml_in_ndm_module, only: rangml
  implicit none

  integer, intent(in)  :: n
  double precision     :: res
  integer  :: i

  if ((n < 0) .and. (rangml == 0)) then
    stop "Argument of factorial function should be positive"
  elseif (n == 0) then
    res = 1
  end if
  res = 1.d0
  do i = 2, n
    res = res*i
  end do
end function factorial



subroutine compute_cg_vector(j_max, ntype)

  use ml_in_ndm_module, only: rangml, cg_vector
  implicit none

  double precision, intent(in)     :: j_max
  integer, intent(in)  :: ntype
  double precision     :: cg
  integer  :: j, j1, j2, m, m1, m2, j0_max, j1_max, j2_max

  if ((ntype == 1) .and. dabs(int(j_max) - j_max) .gt. 0.1d0) then
    if (rangml == 0) then
      write (*, *) 'ML: Warning: probably you intend to use some SO3 or bi-SO3 -  related descriptor'
      write (*, *) 'ML: Warning: in that case l_max can be only integer, now l_max', j_max
    end if
    stop
  end if

  j0_max = int(j_max*ntype)
  j1_max = int(j_max*ntype)
  j2_max = int(j_max*ntype)

  if (allocated(cg_vector)) deallocate (cg_vector)
  ! 6 Dims
  allocate (cg_vector(0:j1_max, -j1_max:j1_max, 0:j2_max, -j2_max:j2_max, 0:j0_max, -j0_max:j0_max))
  cg_vector(:, :, :, :, :, :) = 0.d0
  do j1 = 0, j1_max
    do j2 = 0, j2_max
      do j = 0, j0_max
        do m1 = -j1, j1, ntype
          do m2 = -j2, j2, ntype
            do m = -j, j, ntype
              if (m - m1 - m2 /= 0) cycle
              if (j1 + j2 - j < 0) cycle
              if (j - j1 + j2 < 0) cycle
              if (j - j2 + j1 < 0) cycle
              if (ntype == 2) then
                if (mod(j + j1 + j2, ntype) == 1) cycle
              end if
              cg_vector(j1, m1, j2, m2, j, m) = cg(j1, m1, j2, m2, j, m, ntype)
              ! debug write (12,'(6i5, e20.10)') j1,j2,j,  m1,m2, m, cg(j1,m1,j2,m2,j,m, ntype)
            end do
          end do
        end do
      end do
    end do
  end do
  !stop 'test cg in compute cg vector'
end subroutine compute_cg_vector
