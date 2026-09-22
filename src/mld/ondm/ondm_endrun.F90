
#include "../../MLD_MACROS.INC"

subroutine ondm_arret_ndm()

  use ondm_gen_com_m, only: dmtype
  use mld_logger
  use mld_mpi
  use ml_in_ndm_module, only: ML_MLD_DMTYPE

  implicit none


  if (dmtype /= ML_MLD_DMTYPE) then
    call MPI_FINALIZE(codeml)
  end if

  call log_warning("stop in arret_ndm TODO have to use mld_mpi_finalize or mld_mpi_abort")           ! TODO something else to return in main program
  stop

end subroutine ondm_arret_ndm



subroutine ondm_endrun()

  USE module_kind_variables , ONLY: kind_double
  use ondm_gen_com_m, only: im, lPkbar, iteTemp, itesigma, rang, lprtfat, &
                            angst, it, erg2eV,  dmtype, timel, linstantfda 
  use ondm_tab_imm_m, only: xp, fp  

  ! use posana
  !USE cfg_module
  !use elec_cell, only: sauveelec

  use mld_logger
  use mld_string
  use mld_unit
  use time_measure

  implicit none

  _NAMECURRENT_("ondm_endrun")


  logical  :: ok
  integer  :: i, iunit
  real(kind_double)   :: unitP
  character*5    :: cunitP


  ! if (associated(eatom)) eatom(:)=0  TODO ?

  _MLD_BEGIN_
  if (lPkbar) then
    unitP = 1.0d-9
    cunitP = 'kbar'
  else
    unitP = 1.0
    cunitP = 'd/cm2'
  end if

  ! Un dernier calcul des forces pour la route
  IF (iteTemp .GE. 0) iteTemp = 1
  IF (iteSigma .GE. 0) iteSigma = 1


  ! CALL calfo()

  if (rang == 0) then
    if (lprtfat) then
      iunit = open_file('xifi.dat', 'formatted', action='write')
      write (iunit, '(i0)') im
      write (iunit, '(6g20.8)') (xp(1:3, i)*angst, fp(1:3, i)*erg2eV/angst, i=1, im)
      ok = close_unit(iunit)
    end if
  end if



  if (rang == 0) then

    call log_info('ML: neighbours  time '//vtoa(temps_neigh)//nwl// &
                  'ML: energy      time '//vtoa(temps_energy)//nwl// &
                  'ML: force       time '//vtoa(temps_force)//nwl// &
                  'ML: stress      time '//vtoa(temps_stress)//nwl// &
                  'ML: descriptors time '//vtoa(temps_descripteurs))

    call log_info('--- END OF RUN --- iteration '//vtoa(it)//'time '//vtoa(timel))
  end if

  if (.not. linstantfda) then
    !death_ndm if (iteangle >= 0) call adf()
  end if
  if ((dmtype == 2) .or. (dmtype == 3)) then
    it = 0
  end if

  call ondm_analyse()

  !if ((ldesinteg .EQV. .true.) .and. (itdes == nstepdes)) call desinteg_insert()
  !if (.not. parallele .and. iteanapos >= 0) call anapos(it)

  _MLD_END_
  call ondm_arret_ndm()
  stop

end subroutine ondm_endrun


subroutine ondm_redefine_ty

  use ondm_gen_com_m, only: imm 
  use ondm_var_pot, only: ntyp_buffer, ntyp, ty_buffer, ty, cm_buffer, cm
  use ondm_tab_imm_m, only: ityp, ityp_buffer

  implicit none
  ntyp_buffer = ntyp
  allocate (ityp_buffer(imm), ty_buffer(ntyp), cm_buffer(ntyp))
  ityp_buffer(:) = ityp(:)
  ty_buffer(:) = ty(:)
  cm_buffer(:) = cm(:)

  ntyp = 3
  ityp(1:6) = 2
  ityp(8:15) = 2
  ityp(7) = 3


  deallocate (ty, cm)
  allocate (ty(ntyp), cm(ntyp))

  ty(1) = 'Fe'
  ty(2) = 'Cu'
  ty(3) = 'O '

  cm(1:ntyp) = cm(1)

end subroutine ondm_redefine_ty


subroutine ondm_refix_ty
  use ondm_gen_com_m, only: imm 
  use ondm_var_pot, only: ntyp_buffer, ntyp, ty_buffer, ty, cm_buffer, cm
  use ondm_tab_imm_m, only: ityp, ityp_buffer

  implicit none

  ntyp = ntyp_buffer
  deallocate (ty, ityp)
  allocate (ty(ntyp), cm(ntyp), ityp(imm))
  ityp(:) = ityp_buffer(:)
  ty(:) = ty_buffer(:)
  cm(:) = cm_buffer(:)

end subroutine ondm_refix_ty
