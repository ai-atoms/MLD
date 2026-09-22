
subroutine ondm_caltabi
  ! version du 4 juin 2010, 14h38 - last chaged by MCM
  USE module_kind_variables, ONLY: kind_double
  use ondm_gen_com_m, only: im, imm, rvois, lperiod, lconstrtot,ldemitab, &
                       decal_bc, ldecal_bc, nvperat,nvois, ivoismax, bg, at , &
                       indi, indi2, ncel, last, noxyz, nato 

  use ondm_var_pot, only: ipo , npair, ipotentiel
  use ondm_tab_imm_m, only: xp, ityp, iwmax, iwmax2, ielat
  use ondm_transform_coord, only: ondm_notperiod
  implicit none

  integer  :: iw, iwph, i, ip, j, maxvoi, nvij, iwo
  integer  :: itj, ll
  REAL(kind_double)   :: r2
  real(kind_double), dimension(1:npair) :: rvois2

  real(kind_double), dimension(3) :: xpi, dx, ds

  integer  :: iti, &      ! type de i
              koo, &                  ! cel de i
              ncelvois, ko1, &        ! cel voisine de i
              i1, i2, itemp


  real(kind_double), dimension(:, :), allocatable   :: xpnp                      ! MODIF Cosmin


  ! OUVERTURE BOUCLE SUR I
  iw = 0
  iwph = 0

  if (ipotentiel == 12) then
    rvois2(1) = (7.0d-8)**2
    rvois2(2) = (3.5d-8)**2
    rvois2(3) = (2.8d-8)**2
  else
    rvois2(:) = rvois**2
  end if

  nvij = 0

  ALLOCATE (xpnp(3, imm))
  if (lperiod) then
    xpnp(:, :) = xp(:, :)
  else
    call ondm_notperiod(xp, xpnp)
  end if

  ! construction par kind_double boucle
  if (lconstrtot) then

    do i = 1, im
      iwo = iw
      xpi(:) = xpnp(:, i)
      iti = ityp(i)

      if (ldemitab) then
        ip = i + 1
      else
        ip = 1
      end if

      do j = ip, im
        if (i .eq. j) cycle
        dx(:) = xpi(:) - xpnp(:, j)
        ds(:) = MatMul(dx(:), bg(:, :))
        IF (ldecal_bc .EQV. .TRUE.) THEN
          IF (ds(1) > 0.5) THEN
            ds(3) = ds(3) - decal_bc
          ELSE IF (ds(1) < -0.5) THEN
            ds(3) = ds(3) + decal_bc
          END IF
        END IF                  ! *!

        if ((ds(1) > 0.5d0) .or. (ds(1) < -0.5d0)) ds(1) = ds(1) - dble(Nint(ds(1)))
        if ((ds(2) > 0.5d0) .or. (ds(2) < -0.5d0)) ds(2) = ds(2) - dble(Nint(ds(2)))
        if ((ds(3) > 0.5d0) .or. (ds(3) < -0.5d0)) ds(3) = ds(3) - dble(Nint(ds(3)))

        dx(:) = MatMul(at(:, :), ds(:))
        r2 = Sum(dx(:)**2)

        itj = ityp(j)
        ll = ipo(iti, itj)

        if (r2 > rvois2(ll)) cycle

        iw = iw + 1
        IF (iw .GT. nVois) THEN
          write (0, '(a,i0)') 'Indice iw du tableau de voisin plus grand&
               & que le max, nVois = ', nVois
          write (0, '(a)') 'Augmentez le nombre moyen de voisins par&
               & atome dans le fichier *.din'
          write (0, '(a,i0)') 'valeur actuelle: nvperat = ', nvperat
          STOP '< Caltabi >'
        END IF

        indi(iw) = j
        indi2(iw) = j
      end do
      iwmax(i) = iw
      iwmax2(i) = iw
      nvij = iw - iwo

    end do                  ! im
    maxvoi = iw
    !*************construction par celulle ****************
  else
    !write (*,*) 'THE fist passage'
    do i = 1, im
      iwo = iw
      koo = ielat(i)          ! Numero de la cellule
      iti = ityp(i)
      xpi(:) = xpnp(:, i)

      ncelvois = min(noxyz, 27) - 1
      ! pour chaque cel. voisine
      do i1 = 0, ncelvois
        ko1 = ncel(koo, i1)
        !               write (6,*)'i1 ko1 ',i1,ko1
        if (ko1 == 0) cycle
        loop_j: do i2 = 1, nato(ko1)
          j = last(i2, ko1)
          !                 write (6,*)'j ',j

          if (ldemitab) then
            if (j .le. i) cycle     ! terme deja calcule
          else
            if (j .eq. i) then
              iwph = iwph + 1
              indi2(iwph) = j
              cycle
            end if
          end if

          itj = ityp(j); ll = ipo(iti, itj)

          dx(:) = xpi(:) - xpnp(:, j)
          ds(:) = MatMul(dx(:), bg(:, :))

          IF (ldecal_bc .EQV. .TRUE.) THEN                 ! *!
            IF (ds(1) > 0.5) THEN
              ds(3) = ds(3) - decal_bc
            ELSE IF (ds(1) < -0.5) THEN
              ds(3) = ds(3) + decal_bc
            END IF
          END IF

          if ((ds(1) > 0.5d0) .or. (ds(1) < -0.5d0)) ds(1) = ds(1) - dble(Nint(ds(1)))
          if ((ds(2) > 0.5d0) .or. (ds(2) < -0.5d0)) ds(2) = ds(2) - dble(Nint(ds(2)))
          if ((ds(3) > 0.5d0) .or. (ds(3) < -0.5d0)) ds(3) = ds(3) - dble(Nint(ds(3)))
          dx(:) = MatMul(at(:, :), ds(:))
          r2 = Sum(dx(:)**2)
          !write (*,*) 'ttt', r2, npair 
          if (r2 > rvois2(ll)) cycle
          iw = iw + 1
          iwph = iwph + 1
          !write (6,*)i,koo,ko1,j,iw, at,bg
          indi(iw) = j
          indi2(iwph) = j
        end do loop_j           ! i2

      end do                  ! ncelvois
      iwmax(i) = iw
      iwmax2(i) = iwph
      nvij = iw - iwo
      if (i > 1) then
        itemp = iwmax2(i) - iwmax2(i - 1)
        if (itemp > ivoismax) ivoismax = itemp
      end if


    end do                  ! fin i
    maxvoi = iw


  end if                  ! lconstrtot


  DEALLOCATE (xpnp)

end subroutine ondm_caltabi
