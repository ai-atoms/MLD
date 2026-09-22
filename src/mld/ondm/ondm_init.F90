
#include "../../MLD_MACROS.INC"

subroutine init()
  use ml_in_ndm_module, only: ML_MLD_DMTYPE
  use ondm_gen_com_m, only: tmean, pmean, timel, kinemean, usdh, two,it, dmtype, imd, &
                      im, imf, imana, ltabvois, tstep,  itmax
  use mld_mpi
  use mld_logger
  use mld_string
  use ondm_var_pot, only: npotentiel, ipotentiel, ntrip, typ_and_pot, npair, typ_pot_pair, &
                      lpotentiel, rue_pot, npotmax, rumax, nad, ntyp, na, rue_pot

  implicit none

  _NAMECURRENT_("init")


  integer  :: lufilmpaf, ipotcont
  integer  :: complet = 1 ! flag d'appel a divid : complet : exec de la routine complete

  _MLD_BEGIN_
  tmean = 0.0
  pmean = 0.0
  timel = 0.0
  kinemean = 0.0
  lufilmpaf = 79


  ! allocating the types
  ipotentiel = -1
  ntyp=1
  npair=1
  npair = ntyp*(ntyp + 1)/2; ntrip = ntyp*ntyp*(ntyp + 1)/2
  call ondm_alloc_typ
  typ_pot_pair = 0

  allocate (rue_pot(npotmax))
  rue_pot(:) = 0.

  ! setting the potential
  call log_info('ML: oNDM -- in init -- POTENTIALS -- with ipotential  : '//vtoa(ipotentiel))
  do ipotcont = 0, npotmax
    if (lpotentiel(ipotcont) .EQV. .true.) then
      ipotentiel = ipotcont
    else
      cycle
    end if
    if (ipotentiel .lt. 10) then
    else
      select case (ipotentiel)
      case (20)
        call log_info('ML: oNDM -- in init ---- continue -- with ipotential  : '//vtoa(ipotentiel))  
        call log_info('ML: oNDM -- before  md_init_potential_ml -- MILADY potential')
        call md_init_potential_ml
      end select
    end if
  end do

  usdh = 1/(two*tstep)
  it = 0

  if (dmtype .ne. 9) then
    call ondm_config
  end if

  call ondm_divid(complet)

  call ondm_DynamicalAllocationCell

  do ipotcont = 0, npotmax
    if (lpotentiel(ipotcont) .EQV. .true.) then
      ipotentiel = ipotcont
    else
      cycle
    end if
  end do

  if (npotentiel .gt. 1) then
    call log_info('ML: oNDM init -- in init -- cell split with rumax  : '//vtoa(rumax*1d8))
  end if
  call ondm_neigcel

  ! compcr non pris en charge en parallele
  ! end setting the cell division
  imd = im
  nad(:ntyp) = na(:ntyp)
  imf = im
  imana = min(imd, imd)

  !computing the neighbours for the very first time
  call ondm_caltabt

  if (ltabvois) call ondm_caltabi

  ! MILADY
  if (ipotentiel == 20) then
    call log_info('ML: oNDM init -- in init -- configuration MILADY -- ')
    call log_info('ML: oNDM init -- in init --     dmtype MILADY    -- '//vtoa(dmtype))
  end if

  if (itmax == 0) call ondm_arret_ndm
  call log_info('ML: oNDM init -- in init -- finalization of init --')

  _MLD_END_

end subroutine init
