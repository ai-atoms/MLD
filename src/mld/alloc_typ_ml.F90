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

subroutine alloc_typ_ml()

#if(MLD_NDM)
  use var_pot  ,only:npair,ntyp,cm,ipo,ty,ntrip, rumax
#else
  use ondm_var_pot
  use ondm_gen_com_m, only: dmtype
  use ml_in_ndm_module, only: ML_MLD_DMTYPE
#endif
  implicit none
  integer  :: i, j, k

  npair = ntyp*(ntyp + 1)/2
  ntrip = ntyp*ntyp*(ntyp + 1)/2
  rumax = 0.

#if(MLD_NDM)
!  if (allocated(na)) deallocate (na); allocate (na(ntyp))
  if (allocated(ipo)) deallocate (ipo); allocate (ipo(ntyp, ntyp))
#else
  if (associated(na)) deallocate (na); allocate (na(ntyp))
  if (associated(ipo)) deallocate (ipo); allocate (ipo(ntyp, ntyp))

#endif

  ! initialisation de ipo
  k = 0
  do i = 1, ntyp
    do j = i, ntyp
      k = k + 1
      ipo(i, j) = k
      ipo(j, i) = k
    end do
  end do

#if(MLD_NDM)
  if (allocated(cm)) deallocate (cm); allocate (cm(ntyp))
!  if (allocated(catom)) deallocate (catom); allocate (catom(ntyp))
!!$  if (allocated(ty)) deallocate (ty); allocate (ty(ntyp))
!!$  if (allocated(pot)) deallocate (pot); allocate (pot(4, npair, 0:ngrid + 1))
!!$   if (allocated(pot_d)) deallocate(pot_d); allocate(pot_d(4,npair,0:ngrid+1))
!!$
!!$  if (allocated(q)) deallocate (q); allocate (q(ntyp))
!!$  q(:) = 0
!!$  if (allocated(rc)) deallocate (rc); allocate (rc(ntyp))
!!$  rc(1:ntyp) = rclu(1:ntyp)*1.d-8
!!$  if (allocated(zz)) deallocate (zz); allocate (zz(npair))
!!$  if (allocated(lue_paire)) deallocate (lue_paire); allocate (lue_paire(npair))
!!$  if (allocated(lu_roff_pair)) deallocate (lu_roff_pair); allocate (lu_roff_pair(npair))
!!$  lu_roff_pair(:) = .false.
!!$
!!$  if (allocated(rue_pair)) deallocate (rue_pair); allocate (rue_pair(npair))
!!$  rue_pair(:) = 0.
!!$
!!$  if (allocated(lue_typ)) deallocate (lue_typ)
!!$  allocate (lue_typ(npair))
!!$  if (allocated(typ_pot_pair)) deallocate (typ_pot_pair)
!!$  allocate (typ_pot_pair(npair))
!!$  if (allocated(lue_trip)) deallocate (lue_trip)
!!$  allocate (lue_trip(ntrip))
!!$  lue_trip(:) = .false.; lue_paire(:) = .false.; lue_typ(:) = .false.
!!$  if (allocated(ro)) deallocate (ro)
!!$  allocate (ro(npair))
!!$  if (allocated(dip)) deallocate (dip)
!!$  allocate (dip(npair))
!!$  if (allocated(pm)) deallocate (pm)
!!$  allocate (pm(npair))
!!$  if (allocated(roff1)) deallocate (roff1)
!!$  allocate (roff1(npair))
!!$  if (allocated(roff2)) deallocate (roff2)
!!$  allocate (roff2(npair))
!!$  if (allocated(a_factor)) deallocate (a_factor)
!!$  allocate (a_factor(npair))
!!$  if (allocated(r8p)) deallocate (r8p)
!!$  allocate (r8p(npair))
!!$
!!$  ! if(iterdf.ge.0) then
!!$  if (allocated(coord)) deallocate (coord)
!!$  allocate (coord(ntyp, ntyp, nkmax))
!!$  if (allocated(digr)) deallocate (digr)
!!$  allocate (digr(ntyp, ntyp, nkmax))
!!$  digr = 0.d0
!!$  if (allocated(fda)) deallocate (fda)
!!$  allocate (fda(ntyp, ntyp, ntyp, contmax))
!!$
!!$  if (allocated(nad)) deallocate (nad)
!!$  allocate (nad(ntyp))
!!$  if (allocated(nas)) deallocate (nas)
!!$  allocate (nas(ntyp))
!!$  if (allocated(nai)) deallocate (nai)
!!$  allocate (nai(ntyp))
!end if
#else
  !if (associated(cm)) deallocate (cm); allocate (cm(ntyp))
  !if (associated(catom)) deallocate (catom); allocate (catom(ntyp))
  !if (associated(ty)) deallocate (ty); allocate (ty(ntyp))
  !if (associated(pot)) deallocate (pot); allocate (pot(4, npair, 0:ngrid + 1))
  if (associated(cm))        deallocate(cm);        allocate(cm(ntyp));          cm(:) = 0.0d0
  if (associated(catom))     deallocate(catom);     allocate(catom(ntyp));      catom(:) = 0
  if (associated(ty))        deallocate(ty);        allocate(ty(ntyp));         ty(:) = 'aaa'

  ! pot is a huge array: pot(4, npair, 0:ngrid+1) - skip in pure ML training mode
  if (dmtype /= ML_MLD_DMTYPE) then
    if (associated(pot))     deallocate(pot);       allocate(pot(4,npair,0:ngrid+1)); pot(:,:,:) = 0.0d0
  else
    if (associated(pot))     deallocate(pot)
    nullify(pot)
  end if

  if (associated(q)) deallocate (q); allocate (q(ntyp)) ; q(:) = 0
  if (associated(rc)) deallocate (rc); allocate (rc(ntyp)) ; rc(:) = 0.d0 
  if (associated(rclu)) deallocate (rclu); allocate (rclu(ntyp)) ;  rclu(:) = 0.d0 
  rc(1:ntyp) = rclu(1:ntyp)*1.d-8
  if (associated(zz)) deallocate (zz); allocate (zz(npair)) ; zz(:) = 0.d0 
  if (associated(lue_paire)) deallocate (lue_paire); allocate (lue_paire(npair)) ; lue_paire(:) = .false.
  if (associated(lu_roff_pair)) deallocate (lu_roff_pair); allocate (lu_roff_pair(npair)) ; lu_roff_pair(:) = .false.


  if (associated(rue_pair)) deallocate (rue_pair); allocate (rue_pair(npair)) ; rue_pair(:) = 0.
  if (associated(lue_typ)) deallocate (lue_typ) ; allocate (lue_typ(npair)) ; lue_typ(:)= .false. 
  if (associated(typ_pot_pair)) deallocate (typ_pot_pair) ; allocate (typ_pot_pair(npair)) ; typ_pot_pair(:) = 0
  if (associated(lue_trip)) deallocate (lue_trip) ; allocate (lue_trip(ntrip)) ; lue_trip(:) = .false. 
  
  if (associated(ro)) deallocate (ro)
  allocate (ro(npair))
  ro(:) = 0.d0 
  if (associated(dip)) deallocate (dip)
  allocate (dip(npair))
  dip(:) = 0.0d0 
  if (associated(pm)) deallocate (pm)
  allocate (pm(npair))
  pm(:) = 0.d0 
  if (associated(roff1)) deallocate (roff1)
  allocate (roff1(npair))
  roff1(:) = 0.d0 
  if (associated(roff2)) deallocate (roff2)
  allocate (roff2(npair))
  roff2(:) = 0.d0 
  if (associated(a_factor)) deallocate (a_factor)
  allocate (a_factor(npair))
  a_factor(:) = 0.0d0 
  if (associated(r8p)) deallocate (r8p)
  allocate (r8p(npair))
  r8p(:) = 0.d0 

  ! coord, digr, fda are large classical-potential arrays - skip in pure ML training mode
  if (dmtype /= ML_MLD_DMTYPE) then
    if (associated(coord)) deallocate (coord)
    allocate (coord(ntyp, ntyp, nkmax))
    coord(:,:,:) = 0.0d0 
    if (associated(digr)) deallocate (digr)
    allocate (digr(ntyp, ntyp, nkmax))
    digr(:,:,:) = 0.d0
    if (associated(fda)) deallocate (fda)
    allocate (fda(ntyp, ntyp, ntyp, contmax))
    fda(:,:,:,:)=0.0d0 
  else
    if (associated(coord)) deallocate (coord)
    nullify(coord)
    if (associated(digr)) deallocate (digr)
    nullify(digr)
    if (associated(fda)) deallocate (fda)
    nullify(fda)
  end if

  if (associated(nad)) deallocate (nad)
  allocate (nad(ntyp))
  nad(:) = 0 
  if (associated(nas)) deallocate (nas)
  allocate (nas(ntyp))
  nas(:) = 0 
  if (associated(nai)) deallocate (nai)
  allocate (nai(ntyp))
  nai(:) = 0 
  ! end if
#endif
end subroutine alloc_typ_ml
