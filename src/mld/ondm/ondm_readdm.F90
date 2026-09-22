
#include "../../MLD_MACROS.INC"

subroutine ondm_readdm()

  use ondm_gen_com_m
  use ondm_var_pot

  use mld_mpi
  use mld_mpi_io
  use mld_unit
  use mld_string
  use mld_logger


  implicit none

  _NAMECURRENT_("readdm")


  character(len=:), allocatable    :: tmp, fnamdin, fnamdin_save
  integer  :: ludin, lufilm, lufilmpaf, i, itean, ic, iThermo, itecompcr, ipotcont, iunit
  logical  :: lginread, ltriclin, lpcon, lfissure, tpot, &
              lxFrozen, lyFrozen, lzFrozen, lxyFrozen, lxzFrozen, lyzFrozen, lxyzFrozen, &
              ok, dinok

  ! integer :: imFree     ! nb d'atomes libres


  namelist /input/ itab, itetabvois, itetemp, itesigma, itefcc, itedepla, tdepla, lfilm, &
    tempstop, tempstopcel, dmtype, lFire, ttol, tfroi, itecoordo, tstep, itetimestep, tsfact, &
    tinit, tcooling, tfcou, epcou, lcasca, lfissure, itmax, nitmax, itean, itespebcout, &
    itederive, igen, linstantrdf, iterdf, nrdf, nfda, linstantfda, rclu, itesauv, formatsauv, &
    lrestart, lPathFromGin, tgc, ltabvois, rvois, rskin, ltpcel, nox, noy, noz, imm, dfpred, &
    ltranche, rulayer, iterasmol, lpcon, lprtzlm, pext, wbox, wNose, lpcon2, lpconxyz, tbox, &
    iteangle, ipotentiel, lpotentiel, itesauvposition, itesauvforce, lfilmext, tdepla2, &
    lTcon, Text, iteTconst, lTberendsen, lTNose, lTHoover, nHoover, tauTcon, ldecal_bc, ldyn2D, &
    maxorder, lalea, rsep, &
    h0, sigext, lconstrtot, lEev, lPkbar, deltax, lcorrelvp, lvpread, &
    lcalcjq, dilat, lderive, lTandersen, nuandersen, landerscou, Llangevin, gamlg, ilangevin, &
    lcdp, lsigtyp, ljqbh, lEparat, itebdv, itetemp2, itecompcr, iteanapos, ldislo, epcoudis, &
    fdislo, lnemd, fnemd, fpstop, iseed, fsumstop, sigstop, lcontr, lpr, lUcell, ibordcou, iteplz, nplz, ngrid, lperiod, &
    lprteat, lprteattotm, lprtfat, lprtsigat, lsigatcel, itecfg, npath, nebtype, nebrelaxation, maxneb, kspring, deltaRmax, &
    rcangle, rcrdf, deltaestop, nbmoye, lHcyl, fmt_cin, lginread, ltriclin, nvperat, &
    lFrozen, lxFrozen, lyFrozen, lzFrozen, lxyFrozen, lxzFrozen, lyzFrozen, lxyzFrozen, imFree, imFirstFrozen, &
    natperc, iteanaposneb, ntyp, &
    lbulle, ldesinteg, nstepdes, ides, kspr, xpspr, typspr, tempdes, neb_noise, neb_noise_scale, lsuivinonpbc, lposmoy, &
    eatref, lheat, rheat, iteheat, theat, Eheat, HessianOrder, kappa, niteration, lanczos_step, mdcg_noise_scale, &
    mdcg_noise, lforcetabulate, ivisu, ibound, user_strainrate, user_stress_yz, fdbkcoef, decal_bc, &
    tempdeplainit, debyetemp,  lprtpot,  timemax, tpseuils, lrctest, l2T, depmaxts, tsmin, &
    itesauvinter, units_lammps


  _MLD_BEGIN_
  ! set default values for variables in namelist
  fnamdin = fnam//'.din'

  ! --- variables de dynamique
  imm = 0                 ! dimensionnement des tableaux atomiques
  itab = 10               ! period of cell repartition
  itetabvois = 10         ! periode de calcul de la table des voisins
  tempstop = -1.0         ! temperature of run stop
  tempstopcel = -1.0      ! temperature of run stop
  dmtype = 0
  !dmtype = type of calculation : 1 -> MD
  !                               2 -> quench (trempe) or fire quench
  !                               3 -> gradient conjugue sur les coordonnes cartesiennes
  !                              30 -> gradient conjugue sur les coordonnes reduites
  !                               4 -> Velocity Verlet
  !                               5 -> test des forces
  !                               6 -> analyse des positions en fin de cascade
  !                               7 -> calcul des phonons
  !                               8 -> PR
  !                               9 -> NEB
  !                              10 -> PARIN RAHMAN
  !                              11 -> UN SEUL CALCUL DE FORCES
  !                              12 -> ART
  !                              16 -> SUNDAE
  !                              17 -> MAB
  !                              18 -> ML
  lFire = .false.         ! Fire algorithm is used for quenching (cf tr_fire.F90)
  ttol = 0.0              ! max tolerance for temperature in %
  tfroi = -1.0            ! imposed temperature
  tstep = 1.0             ! timestep in 10^-15 sec unit
  itetimestep = -1        ! period of check in timestep
  tsfact = 10.0           ! change in time step factor
  tinit = -1.0            ! initial temperature
  tcooling = -1.0         ! cooling rate
  tfcou = -1.0            ! temperature of the border of the box
  epcou = -1.0            ! width of the border of the box
  lcasca = .FALSE.        ! cascade Y/N
  lfissure = .FALSE.      ! crack Y/N
  itmax = -1              ! maximum number of iterations
  nitmax = -1             ! maximum number of new iterations after restart
  itederive = -1          ! "derive" correction
  igen = -2               ! type de generation :0 a partir de.gin, +1 a partir de .cin; -1 de gin vers cin puis stop +2 modification de cin puis stop
  lrestart = .FALSE.      ! if T : restarting from an interrupt job
  lPathFromGin = .FALSE.  ! if T : read initial path in gin files *.1.gin, *.2.gin, ... (NEB calculaion)
  tgc = 0.0               ! threshold for CG calculation
  ltabvois = .FALSE.      ! methode de la table des voisins
  lconstrtot = .FALSE.    ! construction de la table des voisins T=kind_double boucle F=via cel.
  rvois = 0.0             ! rayon de la table des voisins
  rskin = 0.4
  ltpcel = .FALSE.        ! output of temperature and stress in each cell
  lforcetabulate = .FALSE.                         ! The derivative of the energy is NOT tabulated. TRUE if it is.
  lginread= .false. 

  itesauv = 1000          ! period for saving
  itesauvposition = 0     ! periode pour sauvegarde des positions en binaire
  itesauvforce = 0        ! periode pour sauvegarde des forces en binaire
  itesauvinter = 0
  formatsauv = 3          ! format of saving always triclin 3 copmplete ; 2 positions only
  fmt_cin = 1             ! format des fichiers .cin 0 : initiale, 1 = para
  dfpred = 0.1            ! eguess for GC calculations and quenching
  nox = -1
  noy = -1
  noz = -1
  ltranche = .FALSE.      ! existence d'une trache gelee
  rulayer = 0.0           ! largeur de la tranche gelee par
  ibordcou = 0            ! refroidissement sur 3 bords ou seuleument z
  lpr = .false.           ! parinnelo rahman a contrainte constante
  sigext = 0.0            ! Symetric tensor related to the external stress

  ! --- Modif Emmanuel Clouet
  h0(1:3, 1:3) = 0.d0     ! Vecteurs de base de la boite de reference en A (Parrinello, Rahman)
  lUcell = .false.        ! affiche l'energie potentielle de la boite
  ! (cela suppose que h0 correspond a l'etat de reference,
  ! ie etat pour laquelle la contrainte est nulle)
  ! ---

  lpcon = .FALSE.         ! algorithm a pression constante a la hache
  lpcon2 = .FALSE.        ! amortissement de la deformation de la boite
  lpconxyz = .FALSE.      ! the relaxation are allowed only along the X, Y and Z axis


  pext = 0.0              ! pression  par defaut
  wbox = 0.0              ! masse de la boite pour Parrinello-Rahman (par defaut egale a 0.5*masse totale
  wNose = 0.0             ! masse de la boite pour thermostat de Nose (par defaut egale a wbox)
  tbox = 1000.0           ! "temps" de la boite
  lTcon = .false.         ! algorithme a temperature constante
  lTberendsen = .false.   ! algorithme a temperature constante
  lTNose = .false.        ! algorithme a temperature constante de Nose
  lTHoover = .false.      ! algorithme a temperature constante de Hoover
  nHoover = 1
  tauTcon = 200.0         ! The rescales "time" for the Berendsen algorithm
  Text = -1.
  iteTconst = itetemp
  lalea = .FALSE.         ! structure initiale aleatoire
  rsep = 1.0              ! Distance de separation pour le tirage aleatoire
  ipotentiel = -1         ! definit type potentiel : 0=Born-Mayer-Huggins, 1=Buckingham, 2=watanabe,3=buck8,4=UO2, 5 terme Morse, 6=SW �πｴﾎｵ縺､� la Vashista ; 7 pot paire tabule ; 10 EAM; 12 ZrC JuLi(+Tersoff Doan)  ; 13 Tersoff coupure COS; 14 Tersoff coupure FD ; 15 tersoff coupure SIN (original) ; 11 Ercollesi ;; -10=LAMMPS atom style atomic; -11 LAMMPS atom style charge (changes only simple.potin)
  npotentiel = 1          ! nb de potentiels
  lpotentiel(:) = .false.
  ntyp = -1               ! le nombre de type DOIT etre specifie si le nombre de potentiel est superieur �πｴﾎｵ縺､� 1
  ! --- PME
  maxorder = 10           ! Ordre du developpement maximal de la PME
  lvpread = .true.
  dilat(:) = 0.0
  lderive = .false.
  lTandersen = .false.    ! temperature constante a la Andersen
  nuandersen = 1.0d14     ! frequence de tirage aleatoire des vitesses en Hz (valeur elevee = pour cascades)
  landerscou = .false.    ! Andersen seulement sur les bords
  lLangevin = .false.     ! Langevin MD
  gamlg = 5d12            ! Gamma deLangevin (= 0.005/1d-15 fera vp*0.995 pour tstep=1d-15)
  ilangevin = 1
  iko = -1
  lcdp = .false.          ! algorithme d'accumulation de defauts ponctuels
  ldislo = .false.        ! calcul de dislocation
  epcoudis = 0.0          ! epaisseur de la couche avec ajout de force pour dislo
  fdislo = 0.0            ! force appliquee aux atomes de bords
  lnemd = .false.         ! Kth par la methode NEMD Evans, P7229
  fnemd = 0
  fpstop = -0.05          ! critere de conv. sur la force par atome max  pour les trempes UNITE = EV/ANG
  sigstop = -0.05         ! critere de conv. sur les contraintes par direction UNITE = kbar
  fsumstop = -0.1         ! critere de conv. sur la force sqrt ( sum_f F_i^2 )  pour les trempes UNITE = EV/ANG
  lcontr = .false.        ! dynamique contrainte (routine contrainte)
  iseed = 0               ! si <>0 controle le tirage aleatoire des vitesses

  ! ---variables d'analyse
  itetemp = 20            ! period of temperature calculation
  itesigma = -1           ! period of stress calculation
  itefcc = -1             ! period of fcc structure analysis
  itedepla = -100         ! period of displacement cal.
  tdepla = 1.0            ! threshold for displacement
  tdepla2 = -1.0          ! second seuil pour calcul des atomes deplaces
  lfilm = .FALSE.         ! film making of displaced atoms
  lfilmext = .FALSE.      ! film par iteration des atomes deplaces
  itecoordo = -100        ! period of coordination calculation
  itean = 0               ! general control for analysis
  iterdf = -1             ! period of RDF calc. : -1 never ; 0 : nrdf last iterations; +iterdf every iterdf iterations
  linstantrdf = .FALSE.   ! F: calculates and prints the average of the RDF, T : calculent RDFrdf iterations
  linstantfda = .FALSE.   ! F: calculates and prints the average of the RDF, T : calcuatent RDFrdf iterations
  nrdf = 0
  nfda = 0

  lprtzlm = .FALSE.       ! plot du nombres d'atomes par tranche suivant z
  lEev = .false.
  lPkbar = .false.
  ! definition des rayons de coupure pour le calcul des coordinences autour de chaque type atomique
  deltax = 0.0
  if (associated(rclu)) deallocate (rclu); allocate (rclu(20))
  rclu(:) = 2.0


  iterasmol = -1          ! <0 --> genere aucun fichier positions pour logiciel rasmol
  iteangle = -1           ! pilote creation de fichier positions pour
  ! >=0 debut et fin d'execution

  lufilm = 89
  lufilmpaf = 79

  lcorrelvp = .false.
  lcalcjq = .false.
  lsigtyp = .false.       ! calul et affichage de la contrainte atomique
  lEparat = .false.       ! calul et affichage de l'energie par atom
  itebdv = -1             ! frequence de calcul des bond valence
  iteplz = 0
  nplz = 1000
  itetemp2 = -1           ! frequence d'ecriture de la temperature dans fichier separe
  itecompcr = -1          ! remplacee par iteanapos
  iteanapos = -1          ! frequence de comparaison avec cristal de reference

  lperiod = .true.        ! conditions periodiques
  lprteat = .false.       ! if you want to print the energy on atom
  lprtsigat = .false.     ! calul et affichage de la contrainte sur chaque atome
  lsigatcel = .false.     ! calul et affichage de la contrainte atomique moyenne sur la cellule
  lprteattotm = .false.   ! energie par atome totale (pot+cin) moyenne
  ngrid = 20000           ! taille de la grille des potentiels
  itecfg = -1             ! ecriture de fichiers .cfg pour AtomEye

  ! inNEB
  nebtype = 2             ! drag methos is the default
  nebrelaxation = 2       ! We relax all the atoms if nebrelaxation==1 only
  ! the most "deplaced" atoms
  npath = 15              ! 15 images of the neb is the default
  maxneb = 700            ! the MAX of NEB steps
  kspring = 1.0           ! the default value for the spring
  deltaRmax = 1.d-2
  neb_noise_scale = 0.001 ! this will affect the 4th digit
  mdcg_noise_scale = 0.001                         ! this will affect the 4th digit
  ! x + x*neb_noise_scale*random,
  ! where "random" is a random number between
  ! 0 and 1
  neb_noise = 0           ! 0 without noise, 1 with noise
  mdcg_noise = 0          ! 0 without noise, 1 with noise
  ! lperiod=.false.             ! pas de conditions periodiques
  ! inNEB

  rcangle = 3.0
  rcrdf = 5.0
  deltaestop = 0.0
  nbmoye = 100
  lHcyl = .false.
  ltriclin = .true.
  lprtfat = .false.

  lFrozen = .FALSE.
  lxfrozen = .FALSE.      ! .true.: certains atomes sont bloques (pas de dynamique)
  lyfrozen = .FALSE.
  lzfrozen = .FALSE.
  lxyfrozen = .FALSE.
  lxzfrozen = .FALSE.
  lyzfrozen = .FALSE.
  lxyzfrozen = .FALSE.
  imFree = -1             ! The index from which all the atoms with the index i >  imFree are frozen. or with i <= imFree are free
  imFirstFrozen = 0       ! The index from which all the atoms with the index i <= imFirstFrozen are frozen the index i >  imFirstFrozen are free
  ! imFirstFree can be used in the same time with imFree

  nvperat = -1            ! nb moyen de voisins par atomes
  natperc = -1
  iteanaposneb = 0
  lbulle = .false.
  ldesinteg = .false.     ! calcul du delta F de la desintegration d'un atome
  nstepdes = -1
  ides = 1
  kspr = 10.0
  xpspr(:) = -1000.
  typspr = 0
  tempdes = -1.0
  lsuivinonpbc = .false.  ! enable or disable a copy of non folded positions (by the pbc conditions) in binary form each itetimestep.
  lposmoy = .false.       ! writes the average position and energy of the atoms in a .mol file
  eatref(:) = 0.

  lheat = .false.
  rheat = 0.
  Theat = 0.0
  Eheat = 0.
  ivisu = 1               ! format de sortie dans rasmol.f90 : ivisu=1=.mol, ivisu=2=vsim mal coded, ivisu=3=xred

  ! management of the specific boundary conditions (free or rigid)
  ibound = 0              ! ibound = 0 : no spe BoundC, ibound = 1 : strain controlled BoundC, ibound = 2 :stress controlled BoundC
  user_strainrate = 0.    ! crystal strainrate (ibound=1)
  user_stress_yz = 0.     ! stress on the surface (ibound=2)
  fdbkcoef = 0.           ! feedback coefficient for the correction of applied stress (ibound=3)
  decal_bc = 0.
  ldecal_bc = .false.
  itespebcout = -1        ! on n'ecrit pas de .cfg pour le film
  ldyn2D = .false.        ! par defaut : bords libres selon Y


  ! --- in SUNDAE
  kappa = 1e6
  niteration = 10000
  lanczos_step = 1.0d-3
  ! ---

  tempdeplainit = -1
  debyetemp = -1
  lprtpot = .false.


  timemax = 1d20
  tpseuils(:) = 0         ! 1:Tmin; 2:abs(T') ; 3: abs(T'') ; 1:abs(P); 2:abs(P') ; 3: abs(P'')

  lrctest = .true.
           ! temperature de coupure pour les pertes electroniques


  l2T = .false.           ! 2T model
  depmaxts = 0.02
  tsmin = 2.0

  units_lammps = 'metal'
  call log_info('file din name = '//fnamdin)

  ! old way
  inquire (file=fnamdin, exist=dinok)

  if (dinok) then 
  ludin = open_old_file(fnamdin, action='read')    ! err=456
  read (ludin, nml=input)
  else
    ipotentiel=20
    imm=127
    igen=0
    dmtype=18 !18 ML         
    lpr=.false.
    nebtype=2
    nebrelaxation=2
    maxneb=400
    kspring = 0.8
    deltaRmax=0.2 
    npath=9
    dfpred=0.1
    tinit=0.0
    tstep=5.0
    fpstop=0.005
    fsumstop=-1.0
    deltax=0.001
    HessianOrder=2
    
    itmax=250000
    ngrid=250000
    itetemp=200
    itesauv=2000
    itedepla=-1
    itecoordo=-1
    lEev=.true. 
    lEparat=.true.
    
    lPkbar=.true.
    itesigma=100
    
    iterasmol=1000000
    itecfg=500000
    
    ltriclin=.true.
    ltabvois=.true.
    lconstrtot=.false.
    lperiod=.false.
    itetabvois=10
    rvois=8.5
    !rvois=12.0
  end if   ! din input. 


  tsmin = tsmin*1d-15
  depmaxts = depmaxts*1d-8
  if ((l2T) .and. (dabs(tsmin - 2.d-15).lt.1.d-40)) then
    tsmin = 2.d-16
    ! depmaxts=0.002
  end if

  timemax = timemax*1d-15


  imm_glob = imm

  if (imm <= 0) call mld_mpi_abort("readdm: zero atoms")
  if (lpcon) call mld_mpi_abort("readdm: lpcon do not exists any more (obsolete)"//nwl// &
                                "use lpr to get safer Parinnello Rahman")
  if (.not. lperiod) call log_info('readdm: no periodic conditions')


  tstep = tstep*utemps
  tauTcon = tauTcon*utemps
  epcou = epcou*1D-8
  tdepla = tdepla*1D-8
  tdepla2 = tdepla2*1D-8
  rsep = rsep*1.0d-8
  if (itean .ne. 0) then
    itetemp = itean
    itesigma = itean
    itedepla = itean
    itecoordo = itean
    lfilm = .true.
    iterasmol = itean
    iterdf = itean
    linstantrdf = .true.
    iteangle = itean
    linstantfda = .true.
  end if


  if (lpkbar .EQV. .true.) then
    pext = pext*1.0d9
    sigext = sigext*1.0d9
  end if


  ! input check
  if (lnemd) then
    lcalcjq = .true.
    if (dabs(fnemd - 0.d0).lt.1.d-100) call mld_mpi_abort("readdm: lnemd = True and fnemd = 0")
    if (ljqbh) call mld_mpi_abort("readdm: lnemd = True and ljqbh = True")
  end if

  if (ljqbh) then
    ! open(unit=59,file='temptranche.mol')
    dmtype = 1
    read (94, *) njqbh, epsil, epcoud, ittherm, ntr
    if (parallele) then
      if (njqbh .ne. 4) call mld_mpi_abort("readdm: njqbh have to be = 4 in parallel mode")
    end if
    if (njqbh == 5) read (94, *) kthg
  end if


  if (lrestart) then
    igen = 1
    call log_info('restart from file')
  end if
  if ((igen < -1) .or. igen > 2) then              ! +1 from file -1 generate then stop 0 generate then run
    call mld_mpi_abort("readdm: wrong igen = "//vtoa(igen))
  end if

  if (itab <= 0) call mld_mpi_abort("readdm: wrong itab = "//vtoa(itab))

  if (parallele) then
    if (ltabvois) then
      ltabvois = .false.
      call log_info('readdm: ltabvois = False in parallel mode')
    end if
    select case (dmtype)
    case (2, 4)
    case default
      call mld_mpi_abort("readdm: parallel mode forbidden with dmtype = "//vtoa(dmtype))
    end select
  end if


#if(DECOUP)
  ltabvois = .false.
  rvois = 0.
#endif

  if ((.not. ltabvois) .and. itab /= 1) then
    !     if (rang == 0) then
    !        write (6, *) ' '
    !        write (6, *) 'Modification obligatoire a itab = 1 '
    !        write (6, *) ' '
    !     endif
    itab = 1
  end if

  if ((itmax .GT. 0) .and. (nitmax .GT. 0)) call log_info("NITMAX PASSE DEVANT ITMAX")
  if ((.not. lrestart) .and. (nitmax .GT. 0)) call mld_mpi_abort('readdm: lrestart False and nitmax > 0')
  if (dmtype .EQ. 11) itmax = 1
  if (itmax < 0) call mld_mpi_abort('readdm: itmax < 0 : '//vtoa(itmax))
  if (itedepla >= 1 .and. tdepla < 0.0) call mld_mpi_abort('readdm: tdepla < 0 : '//vtoa(tdepla))

  !     if ((itedepla.ge.1).and.(tdepla2.lt.0.0)) then
  !       write (6,*)'wrong tdepla2 < 0 '
  !       call arret_ndm
  !     endif

  if (tdepla2 > 0.0 .and. tdepla2 < tdepla) call mld_mpi_abort('readdm: requires tdepla2 >= tdepla : '//vtoa(tdepla2)//vtoa(tdepla))
  if ((lfilm .and. lfilmext) .and. (itedepla < 1)) call mld_mpi_abort('readdm: inconsistency itedepla < 1 with (lfilm and lfilmext) : True')
  if ((itedepla < 1) .and. lfilm .and. (.not. lfilmext)) call mld_mpi_abort('readdm: inconsistency itedepla < 1 with (lfilm and not lfilmext) : True')
  if ((itedepla < 1) .and. lfilmext .and. (.not. lfilm)) call mld_mpi_abort('readdm: inconsistency itedepla < 1 with (not lfilm and lfilmext) : True')

  !      if (ttol==0.0 .and. tfroi<=0.0) then
  !         write (6, *) 'contradiction ttol <-> tfroi '
  !         stop
  !      endif

  if (tstep < 1D-20 .or. tstep > 1D-13) call mld_mpi_abort('readdm: bad time step = '//vtoa(utemps))
  if (tfcou*epcou <= 0.0) call mld_mpi_abort('readdm: tfcou*epcou <= 0 : '//vtoa(tfcou)//vtoa(epcou))
  if (ltriclin .EQV. .false.) call log_warning('readdm: ltriclin False obsolete ?')

  if (dmtype == 9) lprteat = .true.
  if (dmtype == 12) lprteat = .true.
  if (dmtype == 16) lprteat = .true.
  if (lposmoy .EQV. .true.) then
    lprteattotm = .true.
    call log_info('readdm: lposmoy True'//nwl//'posmoyx as average positions and ending output with average energies')
  end if
  if (lprteattotm .EQV. .true.) then
    lprteat = .true.
    call log_info('readdm: lprteattotm True'//nwl//'compute atoms average energies and output minus eatref(eV) (default 0)')
  end if
  if ((lprteattotm .EQV. .true.) .and. (parallele .EQV. .true.)) call mld_mpi_abort('readdm: eattotm an parallel mode unexpected')


  deltax = deltax*A2cm
  if (dmtype == 5 .and. deltax .le. 0) call mld_mpi_abort('readdm: dmtype 5 and deltax <= 0')
  if (dmtype == 7 .and. deltax .le. 0) then
    if (.not. lEev) call mld_mpi_abort('readdm: PHONDY: lEev should be set on .true.'//nwl//'fix as change accordingly')
    call mld_mpi_abort('readdm: dmtype 7 requires deltax > 0'//nwl//'fix as change deltax')
  end if

  if (dmtype == 7) then
    ldemitab = .false.
    if (.not. ((HessianOrder .eq. 1) .or. (HessianOrder .eq. 2) .or. (HessianOrder .eq. 4))) then
      call mld_mpi_abort('readdm: PHONDY: HessianOrder can have only the values 1, 2 or 4'//nwl// &
                         'which corresponds to a Hessian on 2,3 or 5 points'//nwl// &
                         'HessianOrder: '//vtoa(HessianOrder))
    end if
  end if

  if (dmtype == 18) then
    ldemitab = .FALSE.
    call log_info('ML: ldemitab False, compute all pairs')
  end if



  iThermo = 0
  if (lTandersen) iThermo = iThermo + 1
  IF (lTBerendsen) iThermo = iThermo + 1
  IF (lTcon) iThermo = iThermo + 1
  IF (lTNose) iThermo = iThermo + 1
  IF (lTHoover) iThermo = iThermo + 1
  IF (iThermo .GE. 1) THEN
    if (text .le. 0) call mld_mpi_abort('readdm: T Const et Text <= 0')
    if (lTcon) call log_info('constant temperature: lTcon')
    if (lTberendsen) call log_info('constant temperature : lTBerendsen')
    if (lTNose) call log_info('constant temperature: lTNose')
    if (lTHoover) call log_info('constant temperature: lTHoover')
    if (lTandersen) call log_info('constant temperature: lTandersen')
    IF (iThermo .GT. 1) call mld_mpi_abort('readdm: regular thermostat Berendsen/Nose/Hoover/Andersen: requires only one')
  END IF

  if (linstantrdf .and. iterdf == 0) call mld_mpi_abort('readdm: linstantrdf True and iterdf 0 forbidden')
  if (formatsauv > 3 .or. formatsauv < 0) call mld_mpi_abort('readdm: bad formatsauv '//vtoa(formatsauv))
  if (ltabvois .and. (dabs(rvois -0.d0).lt.1.d-100)) call mld_mpi_abort('readdm: bad rvois '//vtoa(rvois))
  rvois = rvois*1.0d-8


  if (lTcon .and. dmtype > 1) call mld_mpi_abort('readdm: lTcon True and dmtype > 1 forbidden')
  if (lTcon .and. lcasca) call mld_mpi_abort('readdm: lTcon True and lacasca True forbidden')
  if (lTcon .and. (Text < 0.)) call mld_mpi_abort('readdm: lTcon True and Text < 0 forbidden')
  if ((lTberendsen) .and. ((dmtype .EQ. 2) .OR. (dmtype .EQ. 3) .OR. (dmtype .EQ. 30))) then
    call mld_mpi_abort('readdm: lTberendsen True and bad dmtype: '//vtoa(dmtype))
  end if
  if (ipotentiel == -1) then
    tpot = .false.
    lpt: do i = 1, npotmax
      if (lpotentiel(i) .EQV. .true.) then
        tpot = .true.
        exit lpt
      end if
    end do lpt
    if (.not. tpot) then
      call log_error("problem ipotentiel: "//vtoa(ipotentiel)//nwl// &
                     "lpotentiel:"//nwl//vtoa(lpotentiel))
      call ondm_endrun()
    end if
  end if
  ! if ((all(lpotentiel)==.false.).and.(ipotentiel==-1)) then ? end if

  if (ipotentiel .ge. 0) lpotentiel(ipotentiel) = .true.
  npotentiel = 0
  do ipotcont = 0, npotmax
    if (lpotentiel(ipotcont) .EQV. .true.) npotentiel = npotentiel + 1
  end do
  if ((lpotentiel(0) .EQV. .true.) .and. (npotentiel .gt. 1)) then
    call mld_mpi_abort('readdm: npotentiel > 1 and lpotentiel(0) True')
  end if

  if ((ntyp == -1) .and. (npotentiel .gt. 1)) then
    call mld_mpi_abort('readdm: npotentiel > 1 and ntyp -1')
  end if
  if ((npotentiel .gt. 1) .and. (lpotentiel(10) .eqv. .true.)) then
    call log_warning('npotentiel > 1 and EAM implies EAM TAB only')
  end if

  if ((lpotentiel(12) .EQV. .true.) .and. (ltabvois .EQV. .true.)) ldemitab = .false.

  if (lrestart .and. lcorrelvp) call mld_mpi_abort('readdm: restart and correlation forbidden')

  if (lcalcjq) then

  end if

  if (itetimestep > 0) then
    if ((dmtype .eq. 1) .or. (dmtype .eq. 2) .or. (dmtype .eq. 4)) then
      call log_info('time step changed each '//vtoa(itetimestep)//'steps')
    else
      call mld_mpi_abort('readdm: itetimestep only with dmtype 1/2/4')
    end if
  end if

  if (lLangevin .and. (Text .le. 0.0)) call mld_mpi_abort('readdm: Langevin and Text <= 0 forbidden')
  if (lLangevin) then
    dmtype = 4
  end if

  ! end check
  if ((lpconxyz) .and. (.NOT. lpr)) then
    call mld_mpi_abort('readdm: lpconxyz can be used only with PR dynamics, or with lpr True')
  end if


  if (lpr) then
    if (dmtype == 2) lprtrp = .true.
    dmtype = 8
    itesigma = 1
    if (dabs(pext - 0.).lt.1.d-100) then
      call log_info('sigext = '//vtoa(sigext))
      call log_info('pext = '//vtoa(pext))
      do ic = 1, 3
        sigext(ic, ic) = sigext(ic, ic) + pext
      end do
      call log_info('new sigext = '//vtoa(sigext))
    end if

    !=== Modif Emmanuel Clouet ================
    ! Verifie si un etat de reference a ete donne
    IF (Sum(h0(1:3, 1:3)**2) .GE. 1.d-30) lUcell = .TRUE.
    ! Transformation A => cm pour le repere de reference
    Pext = (sigext(1, 1) + sigext(2, 2) + sigext(3, 3))/3.d0
    !=== Fin des modifications ================

    h0(1:3, 1:3) = 1e-8*h0(1:3, 1:3)
    ihbox0(:, :) = 1.d0     ! all the dimension of the box can change
    if (lpconxyz) then
      ihbox0(:, :) = 0.d0     ! ALL the dimension are blockef except:
      ihbox0(1, 1) = 1.d0     ! X
      ihbox0(2, 2) = 1.d0     ! Y
      ihbox0(3, 3) = 1.d0     ! Z
    end if

  end if

  if (tcooling > 0) lastcool = 0.0

  ! MPI
  if (itesauvinter .gt. 0) then
    if (mod(itesauvinter, itesauv) .ne. 0) then
      write (6, *) 'itesauvinter n est pas un multiple de intesauv : stop'
      stop
    end if
  end if

  tmp = '-- DM CONFIGURATION --'//nwl

  select case (dmtype)
  case (1)
    tmp = tmp//'DYNAMIQUE MOLECULAIRE VERLET STANDARD'//nwl
  case (2)
    tmp = tmp//'TREMPE RAPIDE'//nwl
  case (3)
    tmp = tmp//'GRADIENT CONJUGUE sur les coordonnees CARTESIENNES'//nwl
  case (30)
    tmp = tmp//'GRADIENT CONJUGUE sur les coordonnees REDUITES'//nwl
  case (4)
    tmp = tmp//' DYNAMIQUE MOLECULAIRE VELOCITY VERLET'//nwl
  case (5)
    tmp = tmp//'TEST DES FORCES'//nwl
  case (7)
    tmp = tmp//'CALCUL DES PHONONS A PARTIR DE POSITIONS DE FORCES NULLES'//nwl
  case (8)
    tmp = tmp//'PARRINELLO RAHMAN AUTOCOHERENT'//nwl
  case (6)
    tmp = tmp//'ANALYSE DES POSITIONS EN FIN DE CASCADE'//nwl
  case (9)
    tmp = tmp//'DRAG OR NEB DYNAMICS'//nwl
    itesauvposition = -1
    itesauvforce = -1
    itetemp = -1
    itesigma = -1
  case (10)
    tmp = tmp//'TREMPE FIRE'//nwl//'Attention les masses atomiques sont toutes celle du type 1'//nwl
  case (11)
    tmp = tmp//'UN CALCUL DE FORCES'//nwl
#if(ART)
  case (12)
    tmp = tmp//'ART nouveau by N. MOUSSEAU'//nwl
#endif
#if(SUNDAE)
  case (16)
    tmp = tmp//'NDM + SUNDAE'//nwl
#endif
#if(MAB)
  case (17)
    tmp = tmp//'NDM + MAB'//nwl
#endif

  case (18)
    tmp = tmp//'NDM + ML'//nwl

  case default
    call mld_mpi_abort('readdm: bad compute type dmtype: '//vtoa(dmtype))
  end select

  call log_info(tmp)

  if (ltranche) then
    call log_warning('jelly slice')
    ! rulayer=rulayer*1.0d-8
    lfrozen = .true.
    if (lcdp .EQV. .true.) call mld_mpi_abort('readdm: jelly slice + DP forbidden')
  end if

  tmp = '-- CONSTRAINTS --'//nwl

  if (lTcon) then
    tmp = tmp//'TEMPERATURE CONSTANTE as main Text: '//vtoa(text)//nwl
  end if
  if (lTBerendsen) then
    tmp = tmp//'TEMPERATURE CONSTANTE as Berendsen Text: '//vtoa(text)//nwl
  end if

  select case (igen)
  case (-1)
    tmp = tmp//'crystal generation only'//nwl
  case (0)
    tmp = tmp//'crystal generation and run'//nwl
  case (1)
    tmp = tmp//'run from file .cin'//nwl
  case (2)
    tmp = tmp//'modification of file .cin'//nwl
  case default
    call mld_mpi_abort('readdm: bad igen: '//vtoa(igen))
  end select

  call log_info(tmp)

  if (ltabvois) then
    if (npotentiel .gt. 1) then
      call mld_mpi_abort('readdm: ltabvois with some potentials not implemented')
      ! demi table ou table complete = prise de tete
    end if

    ! if (tempstopcel.gt.0) ltpcel=.true.
    if (ltpcel) call log_info('pas de contrainte par celulles')

    tmp = "-- NEIGHBOURS --"//nwl
    select case (ipotentiel)
    case (:9)
      ldemitab = .TRUE.
      tmp = tmp//'DEMI-TABLE DES VOISINS rvois '//vtoa(rvois)//nwl
    case (11:18)
      ldemitab = .false.
      tmp = tmp//'DEMI-TABLE DES VOISINS rvois '//vtoa(rvois)//nwl

    case (10)               ! Potentiel EAM
      if (dmtype == 7) then
        ldemitab = .FALSE.
        tmp = tmp//'TABLE DES VOISINS COMPLETE rvois '//vtoa(rvois)//nwl
      else
        ldemitab = .TRUE.
        tmp = tmp//'DEMI-TABLE DES VOISINS rvois '//vtoa(rvois)//nwl
      end if
    case (20)
      if (lconstrtot) then
        call mld_mpi_abort('readdm: MiLaDy potentials should have lconstrtot False with ipotentiel 20')
      end if
    end select

    !     if(ipotentiel.le.10) then
    !        ldemitab=.TRUE.
    !        if (rang.eq.0) write (6, *) '    DEMI-TABLE DES VOISINS rvois ',rvois
    !     else
    !        ldemitab=.false.
    !        if (rang.eq.0) write (6,*)'    TABLE DES VOISINS COMPLETE rvois ',rvois
    !     end if

    call log_info(tmp)

  end if

  if (dmtype == 18) then
    ldemitab = .FALSE.
    call log_info('For MiLady TABLE DES VOISINS COMPLETE rvois '//vtoa(rvois))
    if (lconstrtot) call mld_mpi_abort('readdm: MiLaDy potentials should have lconstrtot set to False')
  end if

  if (lforcetabulate) then
    if (ipotentiel /= 10) then
      call mld_mpi_abort('readdm: There is no implementation for lforcetabulate True and ipotentiel '//vtoa(ipotentiel)//nwl// &
                         'Change lforcetabulate to False or ipotential to EAM (10)')
    end if

  end if

  if (itetemp2 == -1) itetemp2 = itetemp
  !NONEED! tmp = '-- ANALYZES --'//nwl
  !NONEED! tmp = tmp//toline('itetemp', vtoa(itetemp))//nwl
  !NONEED! tmp = tmp//toline('itesigma', vtoa(itesigma))//nwl
  !NONEED! if (itecoordo > 0) tmp = tmp//toline('itecoordo', vtoa(itecoordo))//nwl
  !NONEED! call log_debug(tmp)


  !NONEED! tmp = '-- CONTROLS --'//nwl
  !NONEED! tmp = tmp//toline('tstep', vtoa(tstep))//nwl
  !NONEED! tmp = tmp//toline('itmax', vtoa(itmax))//nwl
  !NONEED! tmp = tmp//toline('itab', vtoa(itab))//nwl
  !NONEED! tmp = tmp//toline('itetimestep', vtoa(itetimestep))//nwl
  !NONEED! if (itederive > 0) tmp = tmp//toline('itederive', vtoa(itederive))//nwl
  !NONEED! tmp = tmp//toline('tinit', vtoa(tinit))//nwl
  !NONEED! if (tempdeplainit .GT. 0) then
  !NONEED!   if (debyetemp == -1) call mld_mpi_abort('readdm: no defined debyetemp but tempdeplainit > 0')
  !NONEED! end if
  !NONEED! tmp = tmp//toline('tempdeplainit', vtoa(tempdeplainit))//nwl
  !NONEED! tmp = tmp//toline('debyetemp', vtoa(debyetemp))//nwl
  !NONEED! 
  !NONEED! if (ttol > 0.) tmp = tmp//toline('ttol', vtoa(ttol))//nwl
  !NONEED! if (tfroi > 0.) tmp = tmp//toline('tfroi', vtoa(tfroi))//nwl
  !NONEED! if ((tempstop > 0.) .and. (rang == 0)) tmp = tmp//toline('tempstop', vtoa(tempstop))//nwl
  !NONEED! if ((tfcou > 0.) .and. (rang == 0)) then
  !NONEED!   tmp = tmp//toline('tfcou', vtoa(tfcou))//nwl
  !NONEED!   tmp = tmp//toline('epcou', vtoa(epcou*1D+8))//nwl
  !NONEED! end if
  !NONEED! if (tcooling > 0.) tmp = tmp//toline('tcooling', vtoa(tcooling))//nwl
  !NONEED! 
  !NONEED! call log_info(tmp)



  if ((lprtsigat .eqv. .true.) .or. (lsigatcel .eqv. .true.)) then
    lsigat = .true.
  else
    lsigat = .false.
  end if

  if (lPrtSigat .and. (.not. ltabvois)) then
    write (6, '(a)') rang, 'contrainte atomique programme en table des voisins&
         & avec un potentiel EAM ou un terme a deux corps seulement'
    call ondm_arret_ndm
  end if







  ! Gestion des atomes bloques
  IF (lFrozen .OR. lxyzFrozen) THEN
    lxFrozen = .true.; lyFrozen = .true.; lzFrozen = .true.
  END IF
  IF (lxyFrozen) THEN
    lxFrozen = .true.; lyFrozen = .true.
  END IF
  IF (lxzFrozen) THEN
    lxFrozen = .true.; lzFrozen = .true.
  END IF
  IF (lyzFrozen) THEN
    lyFrozen = .true.; lzFrozen = .true.
  END IF

  IF (lxFrozen .OR. lyFrozen .OR. lzFrozen) THEN

    IF (((imFree .gt. 0) .AND. (parallele)) .OR. ((imFirstFrozen /= 0) .AND. (parallele))) THEN
      write (0, '(a)') 'Initialisation du tableau free(:) pour&
           & determiner les atomes bloques non implementes en&
           & parallele'
      STOP '< ReadDm >'
    end IF

    IF (dmType .EQ. 8) THEN
      IF (RANG == 0) write (0, '(a)') 'Vous ne pouvez pas utiliser&
           & Parrinello-Rahman tout en maintenant fixes certains&
           & atomes'
      STOP '< ReadDm >'
    END IF

    ! if ((imFree == -1) .and. (rulayer == 0.0) .and. (imFirstFrozen == 0)) then
    !   if (rang == 0) write (6, *) 'LFROZEN+IMFREE=-1 et RULAYER=0 et imFirstFrozen==0 == stop'
    !   stop
    ! end if

    ! Le tableau free controle quels atomes participent a l'energie (utilise par JP a priori)
    ! Le tableau frozen controle quelles coordonnees de quels atomes sont libres de relaxer
    !    i.e. quelles forces doivent etre annulees
    Allocate (Free(1:imm))
    Free(:) = .true.
    Allocate (Frozen(1:3, 1:imm))
    Frozen(:, :) = .false.  ! Tout le monde bouge ... ...

    if (imFree .ne. -1) then
      IF (lxFrozen) Frozen(1, 1 + imFree:imm) = .true. ! x of i>imFree is frozen
      IF (lyFrozen) Frozen(2, 1 + imFree:imm) = .true. ! y of i>imFree is frozen
      IF (lzFrozen) Frozen(3, 1 + imFree:imm) = .true. ! z of i>imFree is frozen
      IF ((rang == 0) .and. (imFree .ne. imm)) THEN
        write (6, '(a)') 'Dynamique / relaxation avec des atomes bloques'
        if (imFirstFrozen == 0) then
          write (6, '(a,i0,a)') "Seuls les atomes d'indice inferieur ou egal a ", &
            imFree, " bougent"
        else
          write (6, '(a,i0,a,i0,a)') "Seuls les atomes d'indice inferieur ou egal a ", &
            imFree, " et superieur et egal a ", imFirstFrozen + 1, "bougent"
        end if
      end IF
    end if



    lFrozen = .true.
  end IF                  ! lxFrozen,lyFrozen,lzFrozen

  if (rulayer .gt. 0.0) then
    rulayer = rulayer*1.0d-8
    IF (rang == 0) write (6, *) 'atomes immobiles fixes par rulayer ', rulayer*1d8
  end if

  ! sortie readdm
  if ((iteanapos .eq. -1) .and. (itecompcr .ne. -1)) iteanapos = itecompcr
  !if (deltaestop .ne. 0) write (6, *) '<<arret apres chnegement de Epot moyenne ', deltaestop, nbmoye

  usdh = 1/(two*tstep)

  ! tmp = ""
  ! if (rang==0) then
  tmp = toline('npotentiel', vtoa(npotentiel))//nwl
  if (npotentiel == 1) then
    tmp = tmp//toline('ipotentiel', vtoa(ipotentiel))//nwl
    ! tmp = tmp // toline('lpotentiel', vtoa(lpotentiel)) //nwl
  else
    do ipotcont = 1, npotmax
      if (lpotentiel(ipotcont) .EQV. .true.) tmp = tmp//toline('active potential', vtoa(ipotcont))//nwl
    end do
    !        if (lcasca.eqv..true.) then
    !           if (rang==0) write (6,*)'ATTENTION!!! npotentiel>1 et ziegler surement faux !!!!'
    !           stop
    !        end if
  end if
  ! end if
  tmp = tmp//toline('fmt_cin', vtoa(fmt_cin))//nwl

  call log_info(tmp)


  rheat = rheat*1d-8




  ! save result complete namelist, for debug/information
  if (mld_rank == 0) then
    fnamdin_save = fnamdin//'.save'
    ! write save namelists
    iunit = open_file(fnamdin_save, form="formatted", action="write")
    write (nml=input, unit=iunit, delim='quote')     ! namelist
    ok = close_unit(iunit)
    ! tmp = repeat(' ', 6000)
    ! write (nml=input, unit=tmp, delim='quote')     ! namelist
    ! tmp = reformat_namelist(tmp)
    tmp = file_read_stream(fnamdin_save, 6000)
    ! call log_debug("readdm input namelist result"//nwl//tmp)
  end if
  _MLD_END_
end subroutine ondm_readdm
