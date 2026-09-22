module ondm_transform_coord
  implicit none 

  contains 

  subroutine ondm_cryst_to_cart(nvec, vec, trmat, iflag)
  ! Version du 01 septembre 2000
  ! transforms the atomic positions or the k-point
  ! components from crystallographic to carthesian coordinates ( iflag=1)
  ! and viceversa ( iflag=-1 ).
  ! Output carth. coordinates are stored in the input ('vec') array.
  !
  ! if iflag=1:
  !    trmat = at ,  basis of the real-space latt.
  !                  for atoms   or
  !          = bg ,  basis of the rec.-space latt.
  !                  for k-points
  ! if iflag=-1: the opposite
  !
  ! Compute the carth. coordinates of each vectors
  ! (atomic positions or k-points components)

  USE module_kind_variables, ONLY: kind_double

  implicit none

  integer, intent(in)  :: nvec
  integer, intent(in)  :: iflag
  real(kind_double), intent(inout)      :: vec(3, nvec)
  real(kind_double), intent(in)   :: trmat(3, 3)

  integer  :: nv
  real(kind_double), dimension(3) :: vau

  do nv = 1, nvec
    if (iflag == 1) then
      vau = trmat(:, 1)*vec(1, nv) + trmat(:, 2)*vec(2, nv) + trmat(:, 3)*vec(3, nv)
    else
      vau = trmat(1, :)*vec(1, nv) + trmat(2, :)*vec(2, nv) + trmat(3, :)*vec(3, nv)
    end if
    vec(:, nv) = vau
    !     Technique anti-bug
    !     Translation des atomes en bord de boite vers xp=0.0
    !        if (iflag.eq.-1) then
    !          do kpol =1,3
    !            if (vec(kpol,nv).ge.1.0) vec(kpol,nv)=0.0
    !          enddo
    !        endif
  end do

  end subroutine ondm_cryst_to_cart



  subroutine ondm_period
  ! version du 09 decembre 2003
  ! Cette routine applique les conditions periodiques par decalage ou +/-1 (triclinique).
  ! Le decalage est fait sur xp xp ET ax.
  ! Ce decalage sur ax permet de mesurer correctement
  ! le deplacements des atomes depuis leur position de depart

  USE module_kind_variables, ONLY: kind_double
  use ondm_gen_com_m, only: ldecal_bc, lperiod, imm, zero, decal_bc, low_limit, &
                       at, bg 
  use ondm_tab_imm_m, only: xp, xpp, ax 

  implicit none

  integer  :: i, ic
  real(kind_double)   :: cpp, xpici
  ! integer,save  :: iperiod


  ! if (rang==0) write (6,*) 'PARA-T entree period'
  if (.not. lperiod) return

  ! iperiod=iperiod+1
  ! write (*,*) 'PBC PBC PBC capitala tarii e  period',iperiod

  IF (ldecal_bc .EQV. .FALSE.) THEN

    call ondm_cryst_to_cart(imm, xp, bg, -1)              ! cart vers cryst
    call ondm_cryst_to_cart(imm, xpp, bg, -1)
    call ondm_cryst_to_cart(imm, ax, bg, -1)
    do i = 1, imm
      do ic = 1, 3
        xpici = xp(ic, i)
        if ((xpici < 0.d0) .OR. (xpici >= 1.d0)) then
          if ((xpici > -low_limit) .and. (xpici < 0.d0)) then
            xp(ic, i) = zero
          else
            cpp = Dble(Floor(xp(ic, i)))
            xpp(ic, i) = xpp(ic, i) - cpp
            ax(ic, i) = ax(ic, i) - cpp
            xp(ic, i) = xpici - cpp
          end if
        end if
      end do
    end do
    call ondm_cryst_to_cart(imm, xp, at, 1)               ! cryst vers cart
    call ondm_cryst_to_cart(imm, xpp, at, 1)
    call ondm_cryst_to_cart(imm, ax, at, 1)

  ELSE IF (ldecal_bc .EQV. .TRUE.) THEN

    call ondm_cryst_to_cart(imm, xp, bg, -1)              ! cart vers cryst
    call ondm_cryst_to_cart(imm, xpp, bg, -1)
    do i = 1, imm
      XP(3, i) = XP(3, i) - DECAL_bc*int(XP(1, i))

      XPP(1, i) = XPP(1, i) - int(XP(1, i))
      XP(1, i) = XP(1, i) - int(XP(1, i))

      XPP(2, i) = XPP(2, i) - int(XP(2, i))
      XP(2, i) = XP(2, i) - int(XP(2, i))

      XPP(3, i) = XPP(3, i) - int(XP(3, i))
      XP(3, i) = XP(3, i) - int(XP(3, i))
    end do
    call ondm_cryst_to_cart(imm, xp, at, 1)               ! cryst vers cart
    call ondm_cryst_to_cart(imm, xpp, at, 1)
  END IF

  ! if (rang==0) write (6,*)'PARA-T sortie period'

  end subroutine ondm_period



  subroutine ondm_notperiod(xp, xpnp)
  ! version du 09 decembre 2003
  ! applique les conditions periodiques par decalage +/-1 (triclinique).
  ! Le decalage est fait sur xp xp ET ax.
  ! Ce decalage sur ax permet de mesurer correctement
  ! le deplacements des atomes depuis leur position de depart

  USE module_kind_variables, ONLY: kind_double
  use ondm_gen_com_m, only: im, imm, at, bg, low_limit, zero 
  implicit none

  real(kind_double), intent(in)   :: xp(3, imm)
  real(kind_double), intent(inout)   :: xpnp(3, imm)

  integer  :: i, ic
  real(kind_double)   :: xpici


  ! iperiod=iperiod+1
  ! write (*,*) 'PBC PBC PBC capitala tarii e ....          ondm_notperiod',iperiod
  xpnp(:, :) = xp(:, :)
  call ondm_cryst_to_cart(imm, xpnp, bg, -1)            ! cart vers cryst

  do i = 1, im
    do ic = 1, 3
      xpici = xpnp(ic, i)

      if ((xpici < 0.d0) .OR. (xpici >= 1.d0)) then
        if ((xpici > -low_limit) .and. (xpici < 0.d0)) then
          ! write (*,*) 'low',low_limit,xpici
          xpnp(ic, i) = zero
        else
          xpnp(ic, i) = xpici - Dble(Floor(xpici))
        end if
      end if

    end do
  end do

  call ondm_cryst_to_cart(imm, xpnp, at, 1)             ! cryst vers cart

  end subroutine ondm_notperiod


  subroutine ondm_recips(a1, a2, a3, b1, b2, b3)
  ! generates the reciprocal lattice vectors b1,b2,b3
  ! given the real space vectors a1,a2,a3. The b's are units of 2 pi/a.

  USE module_kind_variables, ONLY: kind_double

  implicit none

  real(kind_double), intent(in)   :: a1(3)
  real(kind_double), intent(in)   :: a2(3)
  real(kind_double), intent(in)   :: a3(3)
  real(kind_double), intent(out)  :: b1(3)
  real(kind_double), intent(out)  :: b2(3)
  real(kind_double), intent(out)  :: b3(3)

  integer  :: iperm, i, j, k, l, ipol
  real(kind_double)   :: den, s


  ! first we compute the denominator
  den = 0
  i = 1
  j = 2
  k = 3
  s = 1.D0
  do iperm = 1, 3
    den = den + s*a1(i)*a2(j)*a3(k)
    l = i
    i = j
    j = k
    k = l
  end do
  i = 2
  j = 1
  k = 3
  s = -s
  do while (s < 0.D0)
    do iperm = 1, 3
      den = den + s*a1(i)*a2(j)*a3(k)
      l = i
      i = j
      j = k
      k = l
    end do
    i = 2
    j = 1
    k = 3
    s = -s
  end do

  ! here we compute the reciprocal vectors
  i = 1
  j = 2
  k = 3
  do ipol = 1, 3
    b1(ipol) = (a2(j)*a3(k) - a2(k)*a3(j))/den
    b2(ipol) = (a3(j)*a1(k) - a3(k)*a1(j))/den
    b3(ipol) = (a1(j)*a2(k) - a1(k)*a2(j))/den
    l = i
    i = j
    j = k
    k = l
  end do
 end subroutine ondm_recips



real(kind(0.0d0)) function ondm_distmin(a, b)

  USE module_kind_variables, ONLY: kind_double

  implicit none

  real(kind_double), intent(in)   :: a(3)
  real(kind_double), intent(in)   :: b(3)

  real(kind_double), dimension(3) :: vp
  real(kind_double)   :: nvp, na, nb


  na = sqrt(a(1)**2 + a(2)**2 + a(3)**2)
  nb = sqrt(b(1)**2 + b(2)**2 + b(3)**2)
  vp(1) = a(2)*b(3) - a(3)*b(2)
  vp(2) = a(3)*b(1) - a(1)*b(3)
  vp(3) = a(1)*b(2) - a(2)*b(1)
  nvp = sqrt(vp(1)**2 + vp(2)**2 + vp(3)**2)
  ondm_distmin = nvp/max(na, nb)
end function ondm_distmin



real(kind(0.0d0)) function ondm_calcvol(a1, a2, a3)

  USE module_kind_variables, ONLY: kind_double

  implicit none

  real(kind_double), intent(in)   :: a1(3)
  real(kind_double), intent(in)   :: a2(3)
  real(kind_double), intent(in)   :: a3(3)

  integer  :: i, j, k, l, s, iperm


  ondm_calcvol = 0.0
  i = 1
  j = 2
  k = 3
  s = 1.D0
  do iperm = 1, 3
    ondm_calcvol = ondm_calcvol + s*a1(i)*a2(j)*a3(k)
    l = i
    i = j
    j = k
    k = l
  end do
  i = 2
  j = 1
  k = 3
  s = -s
  do while (s < 0.D0)
    do iperm = 1, 3
      ondm_calcvol = ondm_calcvol + s*a1(i)*a2(j)*a3(k)
      l = i
      i = j
      j = k
      k = l
    end do
    i = 2
    j = 1
    k = 3
    s = -s
  end do
  ondm_calcvol = dabs(ondm_calcvol)
end function ondm_calcvol


end module ondm_transform_coord
