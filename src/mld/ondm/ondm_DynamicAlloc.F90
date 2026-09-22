subroutine ondm_DynamicalAllocationCell

  use ondm_gen_com_m, only: ncel, nato, last, deltadist, noxyz, zero, & 
                       natperc 
  !use ondm_var_pot
  implicit none

  if (associated(ncel)) deallocate (ncel); allocate (ncel(0:noxyz, 0:26))
  if (associated(nato))  deallocate(nato);   allocate (nato(0:noxyz))
  if (associated(last))  deallocate(last);  allocate (last(natperc, 0:noxyz))
  if (associated(deltadist))  deallocate(deltadist);  allocate (deltadist(3, 0:26, noxyz))

  !NOWC! lTPcel is always false
  !NOWC! if (lTPcel .EQV. .true.) then
  !NOWC!   if (associated(sigc))  deallocate(sigc) ; allocate (sigc(3, 3, noxyz)); sigc(:, :, :noxyz) = 0.
  !NOWC! end if

  !NOWC! iewald is always 0 
  !NOWC! write(*,*) 'iewald', iewald
  !NOWC! if (iewald .ge. 1) then
  !NOWC!   if (associated(tabv3)) deallocate(tabv3) ;  allocate (tabv3(-ncoucx:ncoucx, -ncoucy:ncoucy, -ncoucz:ncoucz))
  !NOWC!   if (associated(tabf3)) deallocate(tabf3) ;  allocate (tabf3(ntyp, -ncoucx:ncoucx, -ncoucy:ncoucy, -ncoucz:ncoucz))
  !NOWC! end if

  ncel(:noxyz, :26) = zero                         ! et petite initialisation

  !NOWC! L2T is alwaus false   
  !NOWC! if (l2T .eqv. .true.) then
  !NOWC!   if (allocated(tempc))  deallocate(tempc) ; allocate (tempc(noxyz))
  !NOWC!   if (allocated(elossCel))  deallocate(elossCel) ;  allocate (elossCel(noxyz))
  !NOWC! end if
  !NOWC!   lsigatcel is always FALSE 
  !NOWC! if (lsigatcel .eqv. .true.) then
  !NOWC!   if (allocated(patcel)) deallocate (patcel) ; allocate (patcel(noxyz))
  !NOWC!   if (allocated(patcelmax)) deallocate (patcelmax) ; allocate (patcelmax(noxyz))
  !NOWC!   if (allocated(sigatcel)) deallocate (sigatcel) ; allocate (sigatcel(3, 3, noxyz))
  !NOWC!   if (allocated(natchk)) deallocate (natchk) ; allocate (natchk(noxyz))
  !NOWC! end if


end subroutine ondm_DynamicalAllocationCell

subroutine Deallocatecel
  use ondm_gen_com_m
  implicit none
  if (associated(ncel)) deallocate (ncel)
  if (associated(nato)) deallocate (nato)
  if (associated(last)) deallocate (last)
  if (associated(deltadist)) deallocate (deltadist)
  if (associated(sigc)) deallocate (sigc)
  if (associated(tabv3)) deallocate (tabv3)
  if (associated(tabf3)) deallocate (tabf3)
  if (allocated(tempc)) deallocate (tempc)
  if (allocated(tempcm)) deallocate (tempcm)
  if (allocated(tm1)) deallocate (tm1)
  if (allocated(celpp)) deallocate (celpp)
  if (allocated(tcp)) deallocate (tcp)
  if (allocated(pmc)) deallocate (pmc)
  if (allocated(lprtcel)) deallocate (lprtcel)
  if (allocated(elossCel)) deallocate (elossCel)
  if (allocated(patcel)) deallocate (patcel)
  if (allocated(patcelmax)) deallocate (patcelmax)
  if (allocated(sigatcel)) deallocate (sigatcel)
  if (allocated(natchk)) deallocate (natchk)

end subroutine Deallocatecel

subroutine DeallocateVeryAll_varpot 
  use ondm_var_pot
  implicit none 

  !var pot 
  if (associated(na)) deallocate(na)
  !!! if (associated(ipo)) deallocate(ipo)
  if (associated(cm)) deallocate(cm)
  if (associated(cm_buffer)) deallocate(cm_buffer)
  if (associated(catom)) deallocate(catom)
  if (associated(q)) deallocate(q)
  if (associated(rc)) deallocate(rc)
  if (associated(rclu)) deallocate(rclu)
  if (associated(ty)) deallocate(ty)
  if (associated(ty_buffer)) deallocate(ty_buffer)
  if (associated(gamlt)) deallocate(gamlt)
  if (associated(typ_and_pot)) deallocate(typ_and_pot)
  if (associated(typ_pot_pair)) deallocate(typ_pot_pair)
  if (associated(lu_roff_pair)) deallocate(lu_roff_pair)
  if (associated(lue_typ)) deallocate(lue_typ)
  if (associated(lue_trip)) deallocate(lue_trip)
  if (associated(lue_paire)) deallocate(lue_paire)
  if (associated(ipo_2_pair_tab)) deallocate(ipo_2_pair_tab)
  if (associated(pot_pair_tab)) deallocate(pot_pair_tab)
  if (associated(pot)) deallocate(pot)
  if (associated(pot_d)) deallocate(pot_d)
  if (associated(potw)) deallocate(potw)
  if (associated(ray)) deallocate(ray)
  if (associated(shel)) deallocate(shel)
  if (associated(bm)) deallocate(bm)
  if (associated(Dmorse)) deallocate(Dmorse)
  if (associated(amorse)) deallocate(amorse)
  if (associated(remorse)) deallocate(remorse)
  if (associated(zz)) deallocate(zz)
  !!! if (associated(rue_pair)) deallocate(rue_pair)
  !!! if (associated(rue_pot)) deallocate(rue_pot)
  if (associated(gz)) deallocate(gz)
  if (associated(fcr)) deallocate(fcr)
  if (associated(Awat)) deallocate(Awat)
  if (associated(lamb)) deallocate(lamb)
  if (associated(cangle)) deallocate(cangle)
  if (associated(C3C)) deallocate(C3C)
  if (associated(gam)) deallocate(gam)
  if (associated(coup3c)) deallocate(coup3c)
  if (associated(coup3c2)) deallocate(coup3c2)
  if (associated(ipo3c)) deallocate(ipo3c)
  if (associated(l3ctyp)) deallocate(l3ctyp)
  if (associated(l3cpair)) deallocate(l3cpair)
  if (associated(bspg)) deallocate(bspg)
  if (associated(cspg)) deallocate(cspg)
  if (associated(dspg)) deallocate(dspg)
  if (associated(bspf)) deallocate(bspf)
  if (associated(cspf)) deallocate(cspf)
  if (associated(dspf)) deallocate(dspf)
  if (associated(ro)) deallocate(ro)
  if (associated(dip)) deallocate(dip)
  if (associated(pm)) deallocate(pm)
  if (associated(roff1)) deallocate(roff1)
  if (associated(roff2)) deallocate(roff2)
  if (associated(a_factor)) deallocate(a_factor)
  if (associated(r8p)) deallocate(r8p)
  if (associated(bspw)) deallocate(bspw)
  if (associated(cspw)) deallocate(cspw)
  if (associated(dspw)) deallocate(dspw)
  if (associated(capHij)) deallocate(capHij)
  if (associated(capDij)) deallocate(capDij)
  if (associated(capWij)) deallocate(capWij)
  if (associated(ietaij)) deallocate(ietaij)
  if (associated(eamrep)) deallocate(eamrep)
  if (associated(eamrho)) deallocate(eamrho)
  if (associated(eamglue)) deallocate(eamglue)
  if (associated(eamrep_d)) deallocate(eamrep_d)
  if (associated(eamrho_d)) deallocate(eamrho_d)
  if (associated(eamglue_d)) deallocate(eamglue_d)
  if (associated(digr)) deallocate(digr)
  if (associated(coord)) deallocate(coord)
  if (associated(nad)) deallocate(nad)
  if (associated(nas)) deallocate(nas)
  if (associated(nai)) deallocate(nai)
  if (associated(fda)) deallocate(fda)
  if (associated(bsmod1)) deallocate(bsmod1)
  if (associated(bsmod2)) deallocate(bsmod2)
  if (associated(bsmod3)) deallocate(bsmod3)
  if (associated(table)) deallocate(table)
  if (associated(iiim)) deallocate(iiim)
  if (associated(ijim)) deallocate(ijim)
  if (associated(ikim)) deallocate(ikim)
  if (associated(fr1)) deallocate(fr1)
  if (associated(fr2)) deallocate(fr2)
  if (associated(fr3)) deallocate(fr3)
  if (associated(de1)) deallocate(de1)
  if (associated(de2)) deallocate(de2)
  if (associated(de3)) deallocate(de3)

end   subroutine DeallocateVeryAll_varpot

subroutine DeallocateVeryAll_gencomm 
  use ondm_gen_com_m
  implicit none 

!gen_m_com pointer ...
if (associated(xpchup)) deallocate(xpchup)
if (associated(vpchup)) deallocate(vpchup)
if (associated(xpchdn)) deallocate(xpchdn)
if (associated(vpchdn)) deallocate(vpchdn)
if (associated(xpchdeb)) deallocate(xpchdeb)
if (associated(vpchdeb)) deallocate(vpchdeb)
if (associated(itichup)) deallocate(itichup)
if (associated(itichdn)) deallocate(itichdn)
if (associated(itichdeb)) deallocate(itichdeb)
if (associated(num_at_globdesup)) deallocate(num_at_globdesup)
if (associated(num_at_globdesdn)) deallocate(num_at_globdesdn)
if (associated(num_at_globdesdeb)) deallocate(num_at_globdesdeb)
if (associated(Free)) deallocate(Free)
if (associated(Frozen)) deallocate(Frozen)
if (associated(ncel)) deallocate(ncel)
if (associated(nato)) deallocate(nato)
if (associated(last)) deallocate(last)
if (associated(deltadist)) deallocate(deltadist)
if (associated(sigc)) deallocate(sigc)
if (associated(sigat)) deallocate(sigat)
if (associated(sigtyp)) deallocate(sigtyp)
if (associated(sigtyp_loc)) deallocate(sigtyp_loc)
if (associated(sigtyptyp)) deallocate(sigtyptyp)
if (associated(sigtyptyp_loc)) deallocate(sigtyptyp_loc)
if (associated(eatom)) deallocate(eatom)
if (associated(eatomtotm)) deallocate(eatomtotm)
if (associated(indi)) deallocate(indi)
if (associated(indi2)) deallocate(indi2)
if (associated(tabv3)) deallocate(tabv3)
if (associated(tabf3)) deallocate(tabf3)
if (associated(voisins)) deallocate(voisins)
if (associated(latdebord)) deallocate(latdebord)
if (associated(cyl)) deallocate(cyl)
if (associated(b2sINF)) deallocate(b2sINF)
if (associated(b2sSUP)) deallocate(b2sSUP)

  
end subroutine DeallocateVeryAll_gencomm  

subroutine DeallocateAll

  use ondm_gen_com_m
  use ondm_var_pot
  implicit none

  if (associated(ncel)) deallocate (ncel)
  if (associated(nato)) deallocate (nato)
  if (associated(last)) deallocate (last)
  if (associated(deltadist)) deallocate (deltadist)
  if (associated(sigc)) deallocate (sigc)
  if (associated(tabv3)) deallocate (tabv3)
  if (associated(tabf3)) deallocate (tabf3)
  if (associated(na)) deallocate (na)
  !herecos! deallocate (ipo)
  if (associated(cm)) deallocate (cm)

  if (associated(cm_buffer))  deallocate (cm_buffer)
  if (associated(catom))  deallocate (catom)
  if (associated(ty))  deallocate (ty)
  if (associated(ty_buffer))  deallocate (ty_buffer)
  if (associated(gamlt))  deallocate (gamlt)
  if (associated(de1))  deallocate (de1)
  if (associated(de2))  deallocate (de2)
  if (associated(de3))  deallocate (de3)
  if (associated(fr1))  deallocate (fr1)
  if (associated(fr2))  deallocate (fr2)
  if (associated(fr3))  deallocate (fr3)
  if (associated(iiim))  deallocate (iiim)
  if (associated(ijim))  deallocate (ijim)
  if (associated(ikim))  deallocate (ikim)
  if (associated(table))  deallocate (table)
  if (associated(bsmod1))  deallocate (bsmod1)
  if (associated(bsmod2))  deallocate (bsmod2)
  if (associated(bsmod3))  deallocate (bsmod3)
  if (associated(fda))  deallocate (fda)
  if (associated(digr))  deallocate (digr)
  if (associated(ietaij))  deallocate (ietaij)
  if (associated(capDij))  deallocate (capDij)
  if (associated(capHij))  deallocate (capHij)
  if (associated(capWij))  deallocate (capWij)
  if (associated(bspw))  deallocate (bspw)
  if (associated(cspw))  deallocate (cspw)
  if (associated(dspw))  deallocate (dspw)
  !if (associated(rho))  deallocate (rho)
  if (associated(dip))  deallocate (dip)
  if (associated(pm))  deallocate (pm)
  if (associated(dspg))  deallocate (dspg)
  if (associated(bspf))  deallocate (bspf)
  if (associated(cspf))  deallocate (cspf)
  if (associated(bspg))  deallocate (bspg)
  if (associated(dspf))  deallocate (dspf)

  if (associated(pot)) deallocate (pot)
  if (associated(q)) deallocate (q)
  if (associated(rc)) deallocate (rc)
  if (associated(lue_trip)) deallocate (lue_trip)
  if (associated(ro)) deallocate (ro)
  if (associated(roff1)) deallocate (roff1)
  if (associated(roff2)) deallocate (roff2)
  if (associated(a_factor)) deallocate (a_factor)
  if (associated(r8p)) deallocate (r8p)
  if (associated(ray)) deallocate (ray)
  if (associated(bm)) deallocate (bm)
  if (associated(shel)) deallocate (shel)
  if (associated(Awat)) deallocate (Awat)
  if (associated(Bwat)) deallocate (Bwat)
  if (associated(qwat)) deallocate (qwat)
  if (associated(rawat)) deallocate (rawat)
  if (associated(potw)) deallocate (potw)
  if (associated(bspw)) deallocate (bspw)
  if (associated(cspw)) deallocate (cspw)
  if (associated(dspw)) deallocate (dspw)
  if (associated(eamrep)) deallocate (eamrep)
  if (associated(eamrep_d)) deallocate (eamrep_d)
  if (associated(eamglue)) deallocate (eamglue)
  if (associated(eamglue_d)) deallocate (eamglue_d)
  if (associated(eamrho)) deallocate (eamrho)
  if (associated(eamrho_d)) deallocate (eamrho_d)
  if (associated(lamb)) deallocate (lamb)
  if (associated(gam)) deallocate (gam)
  if (associated(cangle)) deallocate (cangle)
  if (associated(coup3c)) deallocate (coup3c)
  if (associated(ipo3c)) deallocate (ipo3c)
  if (associated(coup3c2)) deallocate (coup3c2)
  if (associated(l3ctyp)) deallocate (l3ctyp)
  if (associated(l3cpair)) deallocate (l3cpair)
  if (associated(coord)) deallocate (coord)
  if (associated(digr)) deallocate (digr)
  if (associated(fda)) deallocate (fda)
  if (associated(nad)) deallocate (nad)
  if (associated(nas)) deallocate (nas)
  if (associated(nai)) deallocate (nai)


  if (associated(lue_paire)) deallocate (lue_paire)
  if (associated(lue_typ)) deallocate (lue_typ)
  if (associated(lue_trip)) deallocate (lue_trip)
  if (associated(typ_pot_pair)) deallocate (typ_pot_pair)
  if (associated(typ_and_pot)) deallocate (typ_and_pot)
  if (associated(lu_roff_pair)) deallocate (lu_roff_pair)
  if (associated(pot_pair_tab)) deallocate (pot_pair_tab)
  if (associated(ipo_2_pair_tab)) deallocate (ipo_2_pair_tab)


  
  !herecos! if (associated(rue_pair)) deallocate (rue_pair)
  if (associated(rue_pot)) deallocate (rue_pot)
  if (associated(pot)) deallocate (pot)
  if (associated(potw)) deallocate (potw)
  if (associated(pot_d)) deallocate (pot_d)
  if (associated(ray)) deallocate (ray)
  if (associated(ray)) deallocate (ray)
  if (associated(bm)) deallocate (bm)
  if (associated(Dmorse)) deallocate (Dmorse)
  if (associated(amorse)) deallocate (amorse)
  if (associated(remorse)) deallocate (remorse)

  if (allocated(tempc)) deallocate (tempc)
  if (allocated(tempcm)) deallocate (tempcm)


end subroutine DeallocateAll
