subroutine ondm_param_det()

  USE  module_kind_variables , ONLY: kind_double
  use ondm_gen_com_m, only: rang, zl, pi, lopt, zero 
  use ondm_var_pot, only: npotentiel, ipotentiel, npair, typ_pot_pair, rue_pair, alpha, & 
                     iewald, ncouc3, ncoucx, ncoucy, ncoucz, kpmex,  & 
                     kpmey, kpmez, precis, n2max, nvecttot, kpme 

  implicit none

  integer  :: nc1, l
  real(kind_double)   :: zlm, pparam, k00x, k00y, k00z, k001, rue
  real(kind_double)   :: ruex, ruey, ruez, pparax, pparay, pparaz




  if (npotentiel .ne. 1) then
    if ((iewald .gt. 0) .and. (ncouc3 == 0)) then
      if (rang == 0) write (6, *) 'npot>1 + ewald+ncouc3=0 : stop'
      call ondm_endrun
    end if
    do l = 1, npair
      if ((typ_pot_pair(l) .lt. 10) .and. (typ_pot_pair(l) .ne. 2)) then
        if (rue_pair(l) == 0) then
          if (rang == 0) write (6, *) 'npot>1 + pot paire +ruepaire l =0 : stop', l
          call ondm_endrun
        end if
      end if
    end do

  else
    rue = 0
    if ((ipotentiel .lt. 10) .and. (ipotentiel .ne. 2)) then
      do l = 1, npair
        if (typ_pot_pair(l) .lt. 10) then
          rue = max(rue, rue_pair(l))
        end if
      end do



      ! rue=rue*1.d-8
      alpha = alpha*1.d8
      zlm = max(zl(1), zl(2), zl(3))
      k001 = 2.d0*pi/zlm

      k00x = 2.d0*pi/zl(1)
      k00y = 2.d0*pi/zl(2)
      k00z = 2.d0*pi/zl(3)

      if (ncouc3 == 0 .and. ncoucx /= 0 .and. (ncoucy == 0 .or. ncoucz == 0)) then
        write (6, *) rang, 'Parametres ncouc de la sommation d Exald mal definis'
        call ondm_arret_ndm
      end if

      if (ncouc3 == 0 .and. ncoucy /= 0 .and. (ncoucx == 0 .or. ncoucz == 0)) then
        write (6, *) rang, 'Parametres ncouc de la sommation d Exald mal definis'
        call ondm_arret_ndm
      end if

      if (ncouc3 == 0 .and. ncoucz /= 0 .and. (ncoucx == 0 .or. ncoucy == 0)) then
        write (6, *) rang, 'Parametres ncouc de la sommation d Exald mal definis'
        call ondm_arret_ndm
      end if

      if (ncoucx .lt. 0 .or. ncoucy .lt. 0 .or. ncoucz .lt. 0) then
        write (6, *) rang, 'Parametres ncouc de la sommation d Exald mal definis'
        call ondm_arret_ndm
      end if

      if (ncoucx == 0 .and. ncoucy == 0 .and. ncoucz == 0) then
        nc1 = 0
      else
        nc1 = min(ncoucx, ncoucy, ncoucz)
      end if


      if (.not. lopt) then

        if (precis == zero) then

          if (rue == zero .and. alpha == zero .and. ncouc3 == 0 .and. nc1 == 0) then
            precis = 1.0d-3
            pparam = -log(precis)
            rue = 12.0d-8
            alpha = dsqrt(pparam)/rue
            if (iewald /= 0) then
              ncoucx = int(2.d0*pparam/rue/k00x) + 1
              ncoucy = int(2.d0*pparam/rue/k00y) + 1
              ncoucz = int(2.d0*pparam/rue/k00z) + 1
            end if
            go to 1
          end if

          if (rue == zero .and. alpha == zero .and. ncouc3 /= 0) then
            if (iewald /= 0) then
              precis = 1.0d-3
              pparam = -log(precis)
              ruex = 2.d0*pparam/ncouc3/k00x
              ruey = 2.d0*pparam/ncouc3/k00y
              ruez = 2.d0*pparam/ncouc3/k00z
              rue = max(ruex, ruey, ruez)
              alpha = dsqrt(pparam)/rue
              ncoucx = ncouc3
              ncoucy = ncouc3
              ncoucz = ncouc3
            else
              write (6, *) rang, 'iewald=0 et ncouc3>0 : incoherent !!!'
              call ondm_arret_ndm
            end if
            go to 1
          end if

          if (rue == zero .and. alpha == zero .and. ncouc3 == 0 .and. nc1 /= 0) then
            if (iewald /= 0) then
              precis = 1.0d-3
              pparam = -log(precis)
              ruex = 2.d0*pparam/ncoucx/k00x
              ruey = 2.d0*pparam/ncoucy/k00y
              ruez = 2.d0*pparam/ncoucz/k00z
              rue = max(ruex, ruey, ruez)
              alpha = dsqrt(pparam)/rue
            else
              write (6, *) rang, 'iewald=0 et ncouc/=0 : incoherent !!!'
              call ondm_arret_ndm
            end if
            go to 1
          end if

          if (rue == zero .and. alpha /= zero .and. ncouc3 == 0 .and. nc1 == 0) then
            precis = 1.0d-3
            pparam = -log(precis)
            rue = dsqrt(pparam)/alpha
            if (iewald /= 0) then
              ncoucx = int(2.d0*pparam/rue/k00x) + 1
              ncoucy = int(2.d0*pparam/rue/k00y) + 1
              ncoucz = int(2.d0*pparam/rue/k00z) + 1
            end if
            go to 1
          end if

          if (rue /= zero .and. alpha == zero .and. ncouc3 == 0 .and. nc1 == 0) then
            precis = 1.0d-3
            pparam = -log(precis)
            alpha = dsqrt(pparam)/rue
            if (iewald /= 0) then
              ncoucx = int(2.d0*pparam/rue/k00x) + 1
              ncoucy = int(2.d0*pparam/rue/k00y) + 1
              ncoucz = int(2.d0*pparam/rue/k00z) + 1
            end if
            go to 1
          end if

          if (rue /= zero .and. alpha /= zero .and. ncouc3 == 0 .and. nc1 == 0) then
            pparam = (rue*alpha)**2
            if (iewald /= 0) then
              ncoucx = int(2.d0*pparam/rue/k00x) + 1
              ncoucy = int(2.d0*pparam/rue/k00y) + 1
              ncoucz = int(2.d0*pparam/rue/k00z) + 1
            end if

            go to 1
          end if

          if (rue /= zero .and. alpha == zero .and. ncouc3 /= 0) then
            if (iewald /= 0) then
              pparax = rue*ncouc3*k00x/2.d0
              pparay = rue*ncouc3*k00y/2.d0
              pparaz = rue*ncouc3*k00z/2.d0
              pparam = max(pparax, pparay, pparaz)
              alpha = dsqrt(pparam)/rue
              ncoucx = ncouc3
              ncoucy = ncouc3
              ncoucz = ncouc3
            else
              write (6, *) rang, 'iewald=0 et ncouc3>0 : incoherent !!!'
              call ondm_arret_ndm
            end if
            go to 1
          end if

          if (rue /= zero .and. alpha == zero .and. ncouc3 == 0 .and. nc1 /= 0) then
            if (iewald /= 0) then
              pparax = rue*ncoucx*k00x/2.d0
              pparay = rue*ncoucy*k00y/2.d0
              pparaz = rue*ncoucz*k00z/2.d0
              pparam = max(pparax, pparay, pparaz)
              alpha = dsqrt(pparam)/rue
            else
              write (6, *) rang, 'iewald=0 et ncouc3>0 : incoherent !!!'
              call ondm_arret_ndm
            end if
            go to 1
          end if

          if (rue == zero .and. alpha /= zero .and. ncouc3 /= 0) then
            if (iewald /= 0) then
              ruex = ncouc3*k00x/alpha**2/2.d0
              ruey = ncouc3*k00y/alpha**2/2.d0
              ruez = ncouc3*k00z/alpha**2/2.d0
              rue = max(ruex, ruey, ruez)
              ncoucx = ncouc3
              ncoucy = ncouc3
              ncoucz = ncouc3
            else
              write (6, *) rang, 'iewald=0 et ncouc3>0 : incoherent !!!'
              call ondm_arret_ndm
            end if
            go to 1
          end if

          if (rue == zero .and. alpha /= zero .and. ncouc3 == 0 .and. nc1 /= 0) then
            if (iewald /= 0) then
              ruex = ncoucx*k00x/alpha**2/2.d0
              ruey = ncoucy*k00y/alpha**2/2.d0
              ruez = ncoucz*k00z/alpha**2/2.d0
              rue = max(ruex, ruey, ruez)
            else
              write (6, *) rang, 'iewald=0 et ncouc/=0 : incoherent !!!'
              call ondm_arret_ndm
            end if
            go to 1
          end if

          if (rue /= zero .and. alpha /= zero .and. ncouc3 /= 0) then
            ncoucx = ncouc3
            ncoucy = ncouc3
            ncoucz = ncouc3
            if (iewald == 0) then
              write (6, *) rang, 'iewald=0 et ncouc3>0 : incoherent !!!'
              call ondm_arret_ndm
            end if
            go to 1
          end if

          if (rue /= zero .and. alpha /= zero .and. ncouc3 == 0 .and. nc1 /= 0) then
            if (iewald == 0) then
              write (6, *) rang, 'iewald=0 et ncouc/=0 : incoherent !!!'
              call ondm_arret_ndm
            end if
            go to 1
          end if

          ! endif

        else                    ! (precis /= zero)

          if (rue == zero .and. alpha == zero .and. ncouc3 == 0 .and. nc1 == 0) then
            pparam = -log(precis)
            rue = 12.0d-8
            alpha = dsqrt(pparam)/rue
            if (iewald /= 0) then
              ncoucx = int(2.d0*pparam/rue/k00x) + 1
              ncoucy = int(2.d0*pparam/rue/k00y) + 1
              ncoucz = int(2.d0*pparam/rue/k00z) + 1
            end if

            go to 1
          end if

          if (rue == zero .and. alpha == zero .and. ncouc3 /= 0) then
            if (iewald /= 0) then
              pparam = -log(precis)
              ruex = 2.d0*pparam/ncouc3/k00x
              ruey = 2.d0*pparam/ncouc3/k00y
              ruez = 2.d0*pparam/ncouc3/k00z
              rue = max(ruex, ruey, ruez)
              alpha = dsqrt(pparam)/rue
              ncoucx = ncouc3
              ncoucy = ncouc3
              ncoucz = ncouc3
            else
              write (6, *) rang, 'iewald=0 et ncouc3>0 : incoherent !!!'
              call ondm_arret_ndm
            end if
            go to 1
          end if

          if (rue == zero .and. alpha == zero .and. ncouc3 == 0 .and. nc1 /= 0) then
            if (iewald /= 0) then
              pparam = -log(precis)
              ruex = 2.d0*pparam/ncoucx/k00x
              ruey = 2.d0*pparam/ncoucy/k00y
              ruez = 2.d0*pparam/ncoucz/k00z
              rue = max(ruex, ruey, ruez)
              alpha = dsqrt(pparam)/rue
            else
              write (6, *) rang, 'iewald=0 et ncouc/=0 : incoherent !!!'
              call ondm_arret_ndm
            end if
            go to 1
          end if

          if (rue == zero .and. alpha /= zero .and. ncouc3 == 0 .and. nc1 == 0) then
            pparam = -log(precis)
            rue = dsqrt(pparam)/alpha
            if (iewald /= 0) then
              ncoucx = int(2.d0*pparam/rue/k00x) + 1
              ncoucy = int(2.d0*pparam/rue/k00y) + 1
              ncoucz = int(2.d0*pparam/rue/k00z) + 1
            end if

            go to 1
          end if

          if (rue /= zero .and. alpha == zero .and. ncouc3 == 0 .and. nc1 == 0) then
            pparam = -log(precis)
            alpha = dsqrt(pparam)/rue
            if (iewald /= 0) then
              ncoucx = int(2.d0*pparam/rue/k00x) + 1
              ncoucy = int(2.d0*pparam/rue/k00y) + 1
              ncoucz = int(2.d0*pparam/rue/k00z) + 1
            end if

            go to 1
          end if

          if (rue /= zero .and. alpha /= zero .and. ncouc3 == 0 .and. nc1 == 0) then
            if (rang == 0) &
              write (6, *) 'Seuls la precision et rue sont pris en compte'
            pparam = -log(precis)
            alpha = dsqrt(pparam)/rue
            if (iewald /= 0) then
              ncoucx = int(2.d0*pparam/rue/k00x) + 1
              ncoucy = int(2.d0*pparam/rue/k00y) + 1
              ncoucz = int(2.d0*pparam/rue/k00z) + 1
            end if

            go to 1
          end if

          if (rue /= zero .and. alpha == zero .and. ncouc3 /= 0) then
            if (iewald /= 0) then
              if (rang == 0) &
                write (6, *) 'Seuls la precision et rue sont pris en compte'
              pparam = -log(precis)
              alpha = dsqrt(pparam)/rue
              ncoucx = int(2.d0*pparam/rue/k00x) + 1
              ncoucy = int(2.d0*pparam/rue/k00y) + 1
              ncoucz = int(2.d0*pparam/rue/k00z) + 1
            else
              write (6, *) rang, 'iewald=0 et ncouc3>0 : incoherent !!!'
              call ondm_arret_ndm
            end if
            go to 1
          end if

          if (rue /= zero .and. alpha == zero .and. ncouc3 == 0 .and. nc1 /= 0) then
            if (iewald /= 0) then
              if (rang == 0) &
                write (6, *) 'Seuls la precision et rue sont pris en compte'
              pparam = -log(precis)
              alpha = dsqrt(pparam)/rue
              ncoucx = int(2.d0*pparam/rue/k00x) + 1
              ncoucy = int(2.d0*pparam/rue/k00y) + 1
              ncoucz = int(2.d0*pparam/rue/k00z) + 1
            else
              write (6, *) rang, 'iewald=0 et ncouc/=0 : incoherent !!!'
              call ondm_arret_ndm
            end if
            go to 1
          end if

          if (rue == zero .and. alpha /= zero .and. ncouc3 /= 0) then
            if (iewald /= 0) then
              if (rang == 0) &
                write (6, *) 'Seuls la precision et alpha sont pris en compte'
              pparam = -log(precis)
              rue = dsqrt(pparam)/alpha
              ncoucx = int(2.d0*pparam/rue/k00x) + 1
              ncoucy = int(2.d0*pparam/rue/k00y) + 1
              ncoucz = int(2.d0*pparam/rue/k00z) + 1
            else
              write (6, *) rang, 'iewald=0 et ncouc3>0 : incoherent !!!'
              call ondm_arret_ndm
            end if
            go to 1
          end if

          if (rue == zero .and. alpha /= zero .and. ncouc3 == 0 .and. nc1 /= 0) then
            if (iewald /= 0) then
              if (rang == 0) &
                write (6, *) 'Seuls la precision et alpha sont pris en compte'
              pparam = -log(precis)
              rue = dsqrt(pparam)/alpha
              ncoucx = int(2.d0*pparam/rue/k00x) + 1
              ncoucy = int(2.d0*pparam/rue/k00y) + 1
              ncoucz = int(2.d0*pparam/rue/k00z) + 1
            else
              write (6, *) rang, 'iewald=0 et ncouc/=0 : incoherent !!!'
              call ondm_arret_ndm
            end if
            go to 1
          end if

          if (rue /= zero .and. alpha /= zero .and. ncouc3 /= 0) then
            if (iewald /= 0) then
              if (rang == 0) &
                write (6, *) 'Seuls la precision et rue sont pris en compte'
              pparam = -log(precis)
              alpha = dsqrt(pparam)/rue
              ncoucx = int(2.d0*pparam/rue/k00x) + 1
              ncoucy = int(2.d0*pparam/rue/k00y) + 1
              ncoucz = int(2.d0*pparam/rue/k00z) + 1
            else
              write (6, *) rang, 'iewald=0 et ncouc3>0 : incoherent !!!'
              call ondm_arret_ndm
            end if
            go to 1
          end if

          if (rue /= zero .and. alpha /= zero .and. ncouc3 == 0 .and. nc1 /= 0) then
            if (iewald /= 0) then
              if (rang == 0) &
                write (6, *) 'Seuls la precision et rue sont pris en compte'
              pparam = -log(precis)
              alpha = dsqrt(pparam)/rue
              ncoucx = int(2.d0*pparam/rue/k00x) + 1
              ncoucy = int(2.d0*pparam/rue/k00y) + 1
              ncoucz = int(2.d0*pparam/rue/k00z) + 1
            else
              write (6, *) rang, 'iewald=0 et ncouc/=0 : incoherent !!!'
              call ondm_arret_ndm
            end if
            go to 1
          end if

        end if


1       continue

        n2max = ncouc3*ncouc3

        if (rang == 0) write (6, *) 'Paramtres utilises pour le traitement de EWALD :'

        !determination de kpme
        if (iewald == 2) then
          if (kpmex == 0) kpmex = int(2.5*(2*ncoucx + 1))
          if (kpmey == 0) kpmey = int(2.5*(2*ncoucy + 1))
          if (kpmez == 0) kpmez = int(2.5*(2*ncoucz + 1))
          ! kpmex=100 ; kpmey=100 ; kpmez=100
          kpme = max(kpmex, kpmey, kpmez)
          if (rang == 0) write (6, *) 'kpme x,y,x max ', kpmex, kpmey, kpmez, kpme
        end if

        !          if (kpmex <= 2*ncouc3) then
        !           if (rang==0) write (6,*) 'PME : Taille de grille trop faible STOP!'
        !          stop


        if (rang == 0) then
          if (iewald /= 0) then
            write (6, *) 'RUE=', rue, ' ALPHA=', alpha, ' NCOUC=', ncoucx, ncoucy, ncoucz, &
              ' PRECIS =', precis
          else
            write (6, *) 'RUE=', rue, ' ALPHA=', alpha
          end if
        end if                  ! rang = 0

      else                    ! boucle if (.not. lopt)


#ifdef ixia
        if (rang == 0) write (6, *) 'Optimisation de Ewald pour la machine ixia'
        alpha_ixia = 0.47832d8  ! Valeur optimisee pour ixia

        alpha = alpha_ixia

        if (iewald /= 2) then
          if (rang == 0) write (6, *) 'Demande d optimisation de Ewald et iewald/=2 INCOHERENT !!!'
          stop
        end if

        if (precis == zero) then
          precis = 1.0d-3
          pparam = -log(precis)
          rue = dsqrt(pparam)/alpha_ixia
          ncoucx = int(2.d0*pparam/rue/k00x) + 1
          ncoucy = int(2.d0*pparam/rue/k00y) + 1
          ncoucz = int(2.d0*pparam/rue/k00z) + 1

        else                    ! (precis /= zero)
          pparam = -log(precis)
          rue = dsqrt(pparam)/alpha_ixia
          ncoucx = int(2.d0*pparam/rue/k00x) + 1
          ncoucy = int(2.d0*pparam/rue/k00y) + 1
          ncoucz = int(2.d0*pparam/rue/k00z) + 1

        end if

        if (rang == 0) write (6, *) 'Paramtres utiliss pour le traitement de EWALD(PME) :'

        ! determination de kpme
        kpmex = int(2.5*(2*ncoucx + 1))
        kpmey = int(2.5*(2*ncoucy + 1))
        kpmez = int(2.5*(2*ncoucz + 1))
        kpme = max(kpmex, kpmey, kpmez)
        if (rang == 0) then
          write (6, *) 'kpme x,y,x max ', kpmex, kpmey, kpmez, kpme
          write (6, *) 'RUE=', rue, ' ALPHA=', alpha, ' NCOUC=', ncoucx, ncoucy, ncoucz, ' PRECIS =', precis
        end if

#else
        write (6, *) rang, 'Optimisation de Ewald non prevue pour cette machine'
        call ondm_arret_ndm
#endif


      end if                  ! fin de la boucle if (.not. lopt)

      nvecttot = (2*ncoucx + 1)*(2*ncoucy + 1)*(2*ncoucz + 1) - 1
      write (6, *) 'nvecttot', nvecttot

      rue_pair(:) = rue
    end if
  end if

end subroutine ondm_param_det
