
subroutine ondm_caltabt
  ! Version du 01 fevrier 2001

  USE module_kind_variables, ONLY: kind_double
  use ondm_gen_com_m, only: im, imm, nato, last, lperiod, bg, nox, noy, noz, noxyz, natperc
  use ondm_tab_imm_m, only: xp, ielat 
  use ondm_transform_coord, only: ondm_notperiod, ondm_cryst_to_cart
  implicit none

  integer  :: i, kx, ky, kz, koo
  real(kind_double)   :: aux, auy, auz
  real(kind_double), dimension(:, :), allocatable   :: xpnp                      !


  ! Initialisation
  nato(:noxyz) = 0
  last(natperc, :noxyz) = 0

  ! do i = 1, im
  !    if ((it.ge.1000).and.(i.lt.20)) write (6,'(I5,3G15.7)')i, xp(1,i),xp(2,i),xp(3,i)
  ! end do

  ! cas sans cellule
  if (noxyz == 1) then
    nato(1) = im
    do i = 1, im
      ielat(i) = 1
      last(i, 1) = i
    end do
  else
    ALLOCATE (xpnp(3, imm))
    if (lperiod) then
      xpnp(:, :) = xp(:, :)
    else
      call ondm_notperiod(xp, xpnp)
    end if
    ! Initialisations
    nato(0:noxyz) = 0
    last(natperc, :noxyz) = 0

    ! 1. loop: lattice-coordinates of all atoms
    ! debug write (*,*) 'sub caltabt 1',it,xp(1,1)
    call ondm_cryst_to_cart(imm, xpnp, bg, -1)            ! cart vers cryst
    ! debug write (*,*) 'sub caltabt 2',it,xp(1,1)

    ! if (it.gt.1000) write (6,*) 'CALTABT',it
    do i = 1, im
      ! if  ((it.ge.1000).and.(i.lt.20)) write (6,'(I5,3G15.7)') i, xpnp(1,i),xpnp(2,i),xpnp(3,i)
      aux = xpnp(1, i)*nox
      auy = xpnp(2, i)*noy
      auz = xpnp(3, i)*noz
      kx = int(aux)
      ky = int(auy)
      kz = int(auz)
      ! write (*,*) i, nox,noy,noz, kx,ky,kz
      kx = Modulo(kx, nox)
      ky = Modulo(ky, noy)
      kz = Modulo(kz, noz)
      ! if ((it.ge.1000).and.(i.lt.20)) write (6,'(I5,3G15.7)') i, kx,ky,kz

      koo = 1 + kx + nox*(ky + noy*kz)
      IF ((koo .GT. noxyz) .OR. (koo .LT. 0)) THEN
        write (0, '(a,i0,a,3g20.12)') 'Problem with atom ', i, ', x,y,z = ', xp(1:3, i)
        write (0, '(2(a,i0))') ' koo = ', koo, ' - noxyz = ', noxyz
        STOP
      END IF
      ielat(i) = koo
      nato(koo) = nato(koo) + 1
      ! write (*,*) MAXVAL(nato(:)),koo,i

      ! MODIF Clouet
      IF (nato(koo) .GT. natperc) THEN
        write (0, '(a)') 'You need to increase the maximal number of atoms per cell'
        write (0, '(a,i0)') 'current value: natperc=', natperc
        STOP '< Caltabt >'
      END IF
      ! Fin MODIF Clouet

      last(nato(koo), koo) = i
    end do
    ! debug call cryst_to_cart (imm, xpnp, at, 1)  ! cryst vers cart

    ! do i=1,noxyz
    !   write (6,*) i, nato(i)
    ! end do
    DEALLOCATE (xpnp)       ! MODIF CLOUET
  end if

  ! write (6,*)'sortie caltabt'
  ! write (6,*)'maxnato', maxval(nato)

end subroutine ondm_caltabt
