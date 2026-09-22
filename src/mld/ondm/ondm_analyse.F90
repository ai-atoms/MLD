
#include "../../MLD_MACROS.INC"

subroutine ondm_analyse()

  USE module_kind_variables, ONLY: kind_double
  use ondm_gen_com_m
  use ondm_tab_imm_m
  use former_mat_util, only: matinv
  use ondm_var_pot, only: ntyp, na
 ! USE posana, only: npotmax, lpotentiel, na, l3c, iewald, potisrep, potisglue, ntyp 
  use mld_logger
  use mld_string

  implicit none

  _NAMECURRENT_("analyse")


  integer  :: iti, ic, ko, kx, ky, kz
  real(kind_double), dimension(ntyp)    :: temptyp
  real(kind_double)   :: unitE, unitP
  character*5    :: cunitE, cunitP
  real(kind_double), external     :: tempinst
  integer  :: luvisuc = 888, nprt
  character(len=:), allocatable    :: ext0i
  real(kind_double)   :: xb(3), minp, maxp, mint, maxt
  real(kind_double)   :: minpP, maxpP, mintP, maxtP
  real(kind_double), save   :: Cminp, Cmaxp, Cmint, Cmaxt
  real(kind_double), save   :: CminpP, CmaxpP, CmintP, CmaxtP
  real(kind_double), save   :: CminpP2, CmaxpP2, CmintP2, CmaxtP2
  real(kind_double), save   :: timelm1 = 0

  logical  :: ok
  character(len=:), allocatable    :: tmp
  character(len=*), parameter      :: fmt_136 = '(A, 3I4, 3E15.5, 4E15.7, I4)'


  _MLD_BEGIN_
  !if (allocated(aux_title)) deallocate(aux_title)
  !allocate(aux_title(3))
  !if (allocated(aux_int)) deallocate(aux_int)
  !if (allocated(aux_real)) deallocate(aux_real)

  ! plus //vtoa([itetemp,it]))
  if (it == 1) then
    Cminp = 100000
    Cmaxp = -100000
    Cmint = 100000
    Cmaxt = -100000

    CminpP = 100000
    CmaxpP = -100000
    CmintP = 100000
    CmaxtP = -100000

    CminpP2 = 100000
    CmaxpP2 = -100000
    CmintP2 = 100000
    CmaxtP2 = -100000
  end if
  if (lEev) then
    unitE = erg2eV
    cunitE = '  eV'
  else
    unitE = 1.0
    cunitE = ' erg'
  end if
  if (lPkbar) then
    unitP = 1.0d-9
    cunitP = 'kbar'
  else
    unitP = 1.0
    cunitP = 'd/cm2'
  end if


  if (itetemp > 0) then
    if (mod(it, itetemp) == 0) then

      IF (it <= 1) THEN
        tmean = temp
        kinemean = kine
      ELSE
        tmean = (tmean*(it/itetemp - 1) + temp)/(it/itetemp)
        kinemean = (kinemean*(it/itetemp - 1) + kine)/(it/itetemp)
      END IF

      if (rang == 0) then
        tmp = 'temperature instantanee des atomes de type'//nwl
        do iti = 1, ntyp
          if (na(iti) == 0) cycle
          if (mod(it, itetemp2) == 0) then
            tmp = tmp//toline(vtoa(iti), vtoa(temptyp(iti)))
            ! write (6, '(A,I2,A,F12.2)') '*temp instantanee des atomes de type', iti, ' = ', temptyp(iti)
          end if
        end do
        call log_info(tmp)

        if (mod(it, itetemp2) == 0) then
          ok = .true.
          tmp = '--- stress --- '//cunitP//nwl
        else
          ok = .false.
          tmp = ''                ! precaution
        end if



        if (mod(it, itetemp2) == 0) then
          tmp = '--- valeurs moyennes ---'//nwl//toline('*temperature', vtoa(tmean))//nwl
          ! write (6, '(A,G21.12,A)') '*energie cinetique moyenne = ', kinemean*unitE, cunitE
          if (itesigma > 0) then
            if (mod(it, itesigma) == 0) tmp = tmp//toline('*pression', vtoa(pmean*unitP))
          end if
        end if
        call log_info(tmp)

        if (ltpcel) then
          minp = 100000
          maxp = -100000
          mint = 100000
          maxt = -100000
          !           ext0i = str_0_convi(it, 9)
          !           if (it.ge.3) then
          !              open(luvisuc, file=fnam//'.'//ext0i//'.CEL.mol', form='formatted', status='unknown')
          !              write (luvisuc, '(I9,A,I7,A,D15.6)') 2511 , ' IT =', it, ' Time = ', timel
          !              at=at*1.d8
          !              write (luvisuc,'(9F12.6)')at(1,1),at(2,1),at(3,1),at(1,2),at(2,2),at(3,2),at(1,3),at(2,3),at(3,3)
          !              at=at/1.d8
          !           end if

          nprt = 0
          do kx = 0, nox - 1
            do ky = 0, noy - 1
              do kz = 0, noz - 1
                ko = 1 + kx + nox*(ky + noy*kz)
                xb(1) = float(kx)/float(nox)*at(1, 1) + float(ky)/float(noy)*at(1, 2) + float(kz)/float(noz)*at(1, 3)
                xb(2) = float(kx)/float(nox)*at(2, 1) + float(ky)/float(noy)*at(2, 2) + float(kz)/float(noz)*at(2, 3)
                xb(3) = float(kx)/float(nox)*at(3, 1) + float(ky)/float(noy)*at(3, 2) + float(kz)/float(noz)*at(3, 3)
                xb = xb*1d8
                pmc(ko) = 0.0
                if (itesigma .gt. 0) then
                  if (mod(it, itesigma) .eq. 0) then
                    if (ltabvois) then
                      continue
                    else
                      ! write (6,*)'dans la celulle ',ko,' sigma  '
                      do ic = 1, 3
                        pmc(ko) = pmc(ko) + sigc(ic, ic, ko)/3.0
                        ! write (6,'(3g14.5)')sigc(1,ic,ko),sigc(2,ic,ko),sigc(3,ic,ko)
                        ! write (6,'(A,I2,I2,I2,I2,G14.5,G14.5,G14.5)')'CEL-SIG ',ic,kx,ky,kz,sigc(1,ic,ko),sigc(2,ic,ko),sigc(3,ic,ko)
                      end do
                      ! write (6,'(A,I5,3I4,G14.5,A,A,I4)')'CEL-PRESS ', ko,kx,ky,kz,pmc*unitP, '  ',cunitP,nato(ko)
                      celpP(ko) = 1d-14*unitP*(pmc(ko) - celpm1(ko))/(timel - timelm1)
                      ! celpP2(ko)=1d-14*1d-14*unitP*(pmc(ko)-2*celpm1(ko)+celpm2(ko))/(tstep**2)
                      ! celpp2(ko)=1d-14*1d-14*unitP*((pmc(ko)-celpm1(ko))/tstep -(celpm1(ko)-celpm2(ko))/timelm1)/tstep

                      tcp(ko) = 1d-14*(tempc(ko) - tm1(ko))/(timel - timelm1)
                      ! tcp2(ko)=1d-14*1d-14*(tempc(ko)-2*tm1(ko)+tm2(ko))/(tstep**2)
                      ! tcp2(ko)=1d-14*1d-14*((tempc(ko)-tm1(ko))/tstep -(tm1(ko)-tm2(ko))/timelm1)/tstep
                      ! celpm2(ko)=celpm1(ko)
                      celpm1(ko) = pmc(ko)
                      ! tm2(ko)=tm1(ko)
                      tm1(ko) = tempc(ko)
                      lprtcel(ko) = .false.
                      ! if ((kx==15).or.(ky==15).or.(kz==15).or.(kx==5).or.(ky==5).or.(kz==5)) lprt=.true.
                      ! if((kx.le.13).and.(kx.ge.7).and.(ky.le.13).and.(ky.ge.7).and.(kz.le.13).and.(kz.ge.7)) lprt=.true.
                      ! if (it.le.2) lprt=.false.
                      if (tempc(ko) .gt. tpseuils(1)) lprtcel(ko) = .true.
                      if (abs(tcp(ko)) .gt. tpseuils(2)) lprtcel(ko) = .true.
                      !  if(abs(tcp2(ko)).gt.tpseuils(3)) lprtcel(ko)=.true.
                      if (abs(pmc(ko)*unitP) .gt. tpseuils(3)) lprtcel(ko) = .true.
                      if (abs(celpp(ko)) .gt. tpseuils(4)) lprtcel(ko) = .true.
                      ! if(abs(celpp2(ko)).gt.tpseuils(6)) lprtcel(ko)=.true.
                      if (it .le. 2) lprtcel(ko) = .false.
                      ! lprtcel(ko)=.true.
                      write (6, *) 'tempcprt', tempc(ko), abs(tcp(ko))
                      if (lprtcel(ko) .EQV. .true.) nprt = nprt + 1

                      if (it .ge. 3) then
                        minp = min(minp, pmc(ko)*unitP)
                        maxp = max(maxp, pmc(ko)*unitP)
                        mint = min(mint, tempc(ko))
                        maxt = max(maxt, tempc(ko))

                        minpP = min(minpP, celpp(ko))
                        maxpP = max(maxpP, celpp(ko))
                        mintP = min(mintP, tcp(ko))
                        maxtP = max(maxtP, tcP(ko))
                        ! minpP2=min(minpP2,celpP2(ko))
                        ! maxpP2=max(maxpP2,celpP2(ko))
                        ! mintP2=min(mintP2,tcp2(ko))
                        ! maxtP2=max(maxtP2,tcp2(ko))
                      end if
                    end if
                  end if
                end if
              end do
            end do
          end do
          timelm1 = timel

          if ((mod(it, itetemp2) == 0) .and. (nprt .gt. 0)) then
            write (6, *) '------valeurs par cellules-------', it, nprt
            ext0i = str_0_convi(it, 9)
            open (luvisuc, file=fnam//'.'//ext0i//'.CEL.mol', form='formatted', status='unknown')
            write (luvisuc, '(I9,A,I7,A,D15.6)') nprt, ' IT =', it, ' Time = ', timel
            at = at*1.d8
            write (luvisuc, '(9F12.6)') at(1, 1), at(2, 1), at(3, 1), at(1, 2), at(2, 2), at(3, 2), at(1, 3), at(2, 3), at(3, 3)
            at = at/1.d8
            do kx = 0, nox - 1
              do ky = 0, noy - 1
                do kz = 0, noz - 1
                  ko = 1 + kx + nox*(ky + noy*kz)
                  xb(1) = float(kx)/float(nox)*at(1, 1) + float(ky)/float(noy)*at(1, 2) + float(kz)/float(noz)*at(1, 3)
                  xb(2) = float(kx)/float(nox)*at(2, 1) + float(ky)/float(noy)*at(2, 2) + float(kz)/float(noz)*at(2, 3)
                  xb(3) = float(kx)/float(nox)*at(3, 1) + float(ky)/float(noy)*at(3, 2) + float(kz)/float(noz)*at(3, 3)
                  xb = xb*1d8
                  if (lprtcel(ko) .EQV. .true.) write (luvisuc, fmt_136) 'Au', kx, ky, kz, xb(1), xb(2), &
                    xb(3), tempc(ko), tcp(ko), pmc(ko)*unitP, celpp(ko), nato(ko)
                end do
              end do
            end do
            close (luvisuc)
          end if

          if (it .ge. 3) then
            Cminp = min(Cminp, minp)
            Cmaxp = max(Cmaxp, maxP)
            Cmint = min(Cmint, mint)
            Cmaxt = max(Cmaxt, maxT)

            CminpP = min(CminpP, minpP)
            CmaxpP = max(CmaxpP, maxpP)
            CmintP = min(CmintP, mintP)
            CmaxtP = max(CmaxtP, maxTP)

            ! CminpP2=min(CminpP2,minpP2)
            ! CmaxpP2=max(CmaxpP2,maxPP2)
            ! CmintP2=min(CmintP2,mintP2)
            ! CmaxtP2=max(CmaxtP2,maxTP2)

            write (6, '(A,4G20.10)') 'minmaxp', minp, maxp, mint, maxt
            write (6, '(A,4G20.10)') 'Cminmaxp', Cminp, Cmaxp, Cmint, Cmaxt
            write (6, '(A,4G20.10)') 'minmaxpP', minpP, maxpP, mintP, maxtP
            write (6, '(A,4G20.10)') 'CminmaxpP', CminpP, CmaxpP, CmintP, CmaxtP
            ! write (6,'(A,4G20.10)')'minmaxpP2', minpP2,maxpP2,mintP2,maxtP2
            ! write (6,'(A,4G20.10)')'CminmaxpP2', CminpP2,CmaxpP2,CmintP2,CmaxtP2
          end if
        end if

      end if
    end if
  end if




  _MLD_END_
end subroutine ondm_analyse
