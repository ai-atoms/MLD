subroutine ondm_divid(appel)

  USE module_kind_variables, ONLY: kind_double
  use ondm_gen_com_m, only: itab, ltabvois, rvois, rang, at, bg, volu, normat, &
                            nzl, nox, noy, noz, noxy, noxyz, celsize, zl, im, im_glob, &
                            natperc, nvat, lconstrtot, ldemitab, indi, nvois, pi, indi2, & 
                            nvperat 
  use ondm_var_pot, only: rumax, rue_pair, rue_pot, csive, ngrid, r3cm, r3cm2, &
                          lpotentiel
  !use ondm_tab_imm_m
  use ondm_transform_coord, only: ondm_recips, ondm_calcvol, ondm_distmin

  implicit none


  integer, intent(in)  :: appel                    ! 0: appel partiel juste pour calcul nox/y/z
  ! 1: appel de la routine complete

  real(kind_double)   :: celmin, zlmin, zlm2
  integer  :: izonr2, izonr, ic
  real(kind_double)   :: voluperat, rm2, rut

  call ondm_param_det
  rumax = 0.0
  rumax = max(rumax, maxval(rue_pair))
  csive = rumax/float(ngrid)
  ! write (*,*) rumax, rvois

  ! if(l3c) then
  if (dabs(r3cm - 0.d0).lt.1.d-100)  r3cm = 5.0d-8
  itab = 1
  r3cm2 = r3cm**2
  ! endif

  ! if((rang==0).and.(appel==0)) write (6, *) ' rvois  ', rvois
  if (ltabvois) then
    if (rumax > rvois) then
      if ((rang == 0) .and. (appel == 0)) write (6, '(A,2F12.2)') ' rvois trop petit rvois rumax ', rvois*1d8, rumax*1d8
      ! cos debug call arret_ndm
    else
      rumax = rvois
    end if
  end if

  ! if ((rang==0).and.(appel==0)) write (6, *) 'def rayon coupure ok'

  ! calcul de zlmin
  zlmin = ondm_distmin(at(1, 1), at(1, 2))
  zlm2 = ondm_distmin(at(1, 1), at(1, 3))
  zlmin = min(zlmin, zlm2)
  zlm2 = ondm_distmin(at(1, 2), at(1, 3))
  zlmin = min(zlmin, zlm2)
  zlmin = zlmin*2
  rut = rumax
  if (lpotentiel(10) .eqv. .true.) rut = max(rut, 2*rue_pot(10))
  if (lpotentiel(20) .eqv. .true.) rut = max(rut, 2*rue_pot(20))
  ! write (6,*)'BIP',rumax,rut,rue_pot(10)
  ! end if
  if (lpotentiel(11) .eqv. .true.) rut = max(rut, 2*rue_pot(11))
  if (lpotentiel(12) .eqv. .true.) rut = max(rut, 2*rue_pot(12))
  izonr = int(zlmin/rut)

  ! calcul du volume
  ! if ((rang==0).and.(appel==0)) write (6, *) 'avant volu'
  volu = ondm_calcvol(at(1, 1), at(1, 2), at(1, 3))


  !DETERMINATION DE NOX NOY NOZ


  ! MODIF CLOUET 1
  ! nzl(1:3) doivent etre calcules ici: ils etaient calcules apres l'appel a divid dans config, ce qui n'etait pas correct
  call ondm_recips(at(1, 1), at(1, 2), at(1, 3), bg(1, 1), bg(1, 2), bg(1, 3))
  do ic = 1, 3
    normat(ic) = sqrt(sum(bg(:, ic)**2))
    nzl(ic) = 1.0/normat(ic)
  end do
  ! FIN MODIF CLOUET 1

  if (nox <= 0 .or. noy <= 0 .or. noz <= 0) then
    ! determination de nox noy noz qui ne sont pas donnes dans .din

    ! MODIF CLOUET 2
    nox = int(nzl(1)/rumax)
    noy = int(nzl(2)/rumax)
    noz = int(nzl(3)/rumax)

    IF (nox .LT. 3) nox = 3
    IF (noy .LT. 3) noy = 3
    IF (noz .LT. 3) noz = 3

    IF ((nox .LE. 3) .AND. (noy .LE. 3) .AND. (noz .LE. 3)) THEN
      nox = 1; noy = 1; noz = 1
      ! ltabvois=.TRUE.
      ! lconstrtot=.TRUE.
    END IF

    celsize(1) = zl(1)/float(nox)
    celsize(2) = zl(2)/float(noy)
    celsize(3) = zl(3)/float(noz)

  else

    ! nox noy noz sont donnes dans.din
    if (nox == 2 .or. noy == 2 .or. noz == 2) then
      write (6, *) rang, 'wrong noxyz stop'
      call ondm_arret_ndm
    end if


    celsize(1) = zl(1)/float(nox)
    celsize(2) = zl(2)/float(noy)
    celsize(3) = zl(3)/float(noz)
    celmin = min(celsize(1), celsize(2))
    celmin = min(celsize(3), celmin)

  end if

  ! nox noy et noz sont determines
  noxy = nox*noy
  noxyz = nox*noy*noz

  if (appel == 0) then
    return
  end if


  IF (natperc .LE. 0) THEN                         ! MODIF Clouet
    natperc = INT(im_glob/noxyz)
    nvat = 10*natperc
    natperc = max(5*natperc, 10)                     ! MODIF Clouet
  ELSE                    ! MODIF Clouet
    nvat = 10*natperc       ! MODIF Clouet
  END IF                  ! MODIF Clouet



  ! natperc= INT(im_glob/noxyz)
  ! nvat=10*natperc

  natperc = max(5*natperc, 20)

  if (ltabvois) then
    if (rumax > rvois) then
      write (6, *) rang, ' rvois trop petit rvois rumax ', rvois, rumax
      call ondm_arret_ndm
    end if
    ! crc rm2=max(rumax,2*rvois)
    rm2 = max(rumax, rvois)
    izonr2 = int(zlmin/rm2)
    if (izonr2 < 1) then
      write (6, *) rang, izonr2, rumax, rvois, 'trop petite boite pour rvois !!!'
      !cosboite   call arret_ndm
      write (6, *) rang, 'trop petite boite pour rvois !!!'
    end if

    if (.not. lconstrtot) rumax = rvois
    ! write (*,*) 'DEBUG IN DIVID volu, im', volu, im
    voluperat = volu/im
    IF (nvperat .LE. 0) nvperat = int(4.d0*Pi*(rvois + 1.0d-8)**3/(3*voluperat))


    if (ldemitab) then
      nvois = max(Int(1.5*nvperat*im), 100)
      nvat = max(Int(nvperat*1.3), 10)
    else
      nvois = max(Int(1.5*nvperat*im), 100)
      nvat = max(Int(nvperat*1.3), 10)
    end if


    allocate (indi(nvois))
    allocate (indi2(nvois))
  end if !ltabvois 

  izonr = int(zlmin/r3cm)

end subroutine ondm_divid
