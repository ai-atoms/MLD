import globalv
import os

def init():
 global EXELAMMPS, ELASTIC, eampot, lammps_exe 
 global dinfile, dinfile_neb,  dinfile_ph
 global base, inpdir
 global nodes, ntasks, type_node, nkpar, dynamics
 #potentially used in furture .... 
 global dinfile_ph, dinfile_cg, dinfile_cg2,dinfile_tr
 global fblock, ffreq, fu, fv, fm
 global flmodes, fldos,eigenvalues_phond
 global input_y
 nrepbox=4.0 
 base=os.getcwd()
 inpdir=base+'/INP_NDM/' 
 prefix=globalv.prefix
 nodes=globalv.nodes
 ntasks=globalv.ntasks 
 nkpar=globalv.npools
 type_node=globalv.type_node
 dynamics=globalv.dynamics
 lammps_exe=globalv.lammps_exe


 # occigen, zelda, my_computer
 computer=globalv.cluster

 if computer=='occigen':
   ELASTIC=os.environ['HOME']+'/Elastic/'
   EXELAMMPS=os.environ['STOREDIR']+'/wPH/NDM/bin/rundm90'
   eampot=os.environ['STOREDIR']+"/wPH/eamtab.potin_M2EX_25p_new"
 if computer=='zelda':
   ELASTIC=os.environ['HOME']+'/FitFormation/Elastic/'
   #EXELAMMPS=os.environ['HOME']+'/NDM/bin/rundm90'
   EXELAMMPS='/soft/Lammps/Bin/lmp_mpi'
   #eampot=os.environ['HOME']+"/Flow/INPUT/eamtab.potin_M2EX_25p_new"
   #eampot=os.environ['HOME']+'/Flow/INPUT/Fe_AM05.fs'
   eampot=os.environ['HOME']+'/Flow/INPUT/Al_Liu_2004.fs'
 if computer=='my_computer':
   ELASTIC=os.environ['HOME']+'Elastic/'
   #EXELAMMPS=os.environ['HOME']+'/Codes/lammps-31Mar17/src/lmp_serial'
   #EXELAMMPS=os.environ['HOME']+'/Flow/lammps_milady.git/src/lmp_mpi'
   #EXELAMMPS='/ccc/work/cont002/den/zhonganr/Flow/milady_lammps_beta/src/lmp_mpi'
   #EXELAMMPS=os.environ['HOME']+'/GitHub/mcm_milady_lammps.git/src/lmp_mpi '
   #EXELAMMPS=os.environ['CCCWORKDIR']+'/Flow/mcm_milady_lammps.git/src/lmp_mpi '
   eampot=os.environ['HOME']+"/Codes/lammps-31Mar17/potentials/Al_Liu_2004.fs"
   
 EXELAMMPS=lammps_exe 
 
 if not( (dynamics=='min') or ( dynamics=='neb')):
   print('dynamics can be only min or polar. Now is: %s'%dynamics)
   exit(0)
 
 root_dir=globalv.root_dir
 dirNDM=root_dir+'/'+globalv.input_dirLAMMPS
 dinfile=dirNDM+'/'+globalv.input_incar
 eampot=globalv.input_potcar
 if dynamics=='neb':
   dinfile_neb=dirNDM+'/'+globalv.input_neb_incar
   if not os.path.exists(dinfile_neb):
     print("lammps setup_lammps: in file %s doesn't exist. "%str(dinfile_neb))
     exit(0)

   if ntasks  != 12:
     print("lammps setup_lammps: the number of procs should be  exactely 12 for neb calculations")
     exit(0)

 # relaxation = 'yes' / 'no'
 # if the relaxation is performed.E.g. no - in the case of bulk structures; yes - for SIAs or others unrelaxed strunctures
 

 if not os.path.exists(dinfile):
    print("lammps setup_lammps: in file %s doesn't exist. Put the correct path in setup_lammps"%str(dinfile))
    exit(0)

 if not os.path.exists(eampot):
    print("lammps setup_lammps: Potential file %s doesn't exist. Put the correct path in setup_lammps"%str(eampot))
    exit(0)

 script=globalv.script
 cal_elas=globalv.cal_elas
 if cal_elas=='qha_bulk':
    print("not yet implementation for qha and lammps mode.")
    exit(0)
  
    #limit for zero in frequencies (THz untis).  In order to detect  the three PBC modes
    limit_omega=globalv.limit_omega
    relaxation=globalv.relaxation
    #should be nodes multiply by 24 ntasks-per-node
    nodes=globalv.nodes
    ntasks=globalv.ntasks
    eigenvalues_phondy="eigenvalues.dat"
    dinfile_ph=dirNDM+'/PHONDY/cube.din_ph'
    if computer=='zelda':
      EXENDM_PH=os.environ['HOME']+ '/NDM/p2.8/rundm90_phondy_para'
      SAUVE2GIN=os.environ['HOME']+ '/bin/sauve2gin'
    if computer=='my_computer':
      EXENDM_PH=os.environ['HOME']+ '/NDM/p2.8/rundm90_phondy_para_lammps'
      SAUVE2GIN=os.environ['HOME']+ '/bin/sauve2gin'
    if computer=='occigen':
      EXENDM_PH=os.environ['STOREDIR']+ '/wPH/NDM/bin/rundm90_phondy_para'
      SAUVE2GIN=os.environ['HOME']+ '/bin/sauve2gin'
    if not os.path.exists(dinfile_ph):
       print("ndm setup_ndm_ph: DIN PH file %s doesn't exist. \
          Put the correct path in setup_ndm_ph"%str(dinfile_ph))
       exit(0)
    flmodes=dirNDM+'/PHONDY/cube.phondy.lmodes'
    fldos  =dirNDM+'/PHONDY/cube.phondy.ldos'
    inpmab =dirNDM+'/PHONDY/cube.phondy'

    dinfile_cg=dirNDM+'/PHONDY/cube.din_cg'
    dinfile_cg2=dirNDM+'/PHONDY/cube.din_cg2'
    dinfile_tr=dirNDM+'/PHONDY/cube.din_tr'


    for testf in (EXENDM_PH, dinfile_ph, inpmab, flmodes, fldos, dinfile_cg, dinfile_cg2,dinfile_tr):
      if not os.path.exists(testf):
        print("ndm setup_ndm_ph: The file %s doesn't exist. Put the correct path. "%(str(testf)))
        exit(0)

 return ;
