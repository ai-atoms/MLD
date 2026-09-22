import globalv
import os

def init():
 global EXENDM, EXENDM_PH, EXEMIN_PH, ELASTIC, SAUVE2GIN, eampot
 global dinfile, mlfile, dinfile_ph
 global dinfile_cg, dinfile_cg2,dinfile_tr
 global base, inpdir
 global mdsub_o, inpmab, namein, ginin
 global fblock, ffreq, fu, fv, fm
 global flmodes, fldos,eigenvalues_phondy
 global nodes, ntasks
 global minexe_ph, mininput_ph, minpot01_ph, minpot02_ph
 nrepbox=4.0
 base=os.getcwd()
 inpdir=base+'/INP_NDM/'

 # occigen, zelda, my_computer
 computer=globalv.cluster

 if computer=='occigen':
   ELASTIC=os.environ['HOME']+'/Elastic.git/'
   #EXENDM=os.environ['HOME']+'NDM/bin/rundm90_ml_para'
   #eampot=os.environ['STOREDIR']+"/wPH/eamtab.potin_M2EX_25p_new"
   EXENDM=globalv.ndm_exe
 if computer=='zelda':
   ELASTIC=os.environ['HOME']+'/FitFormation/Elastic/'
   EXENDM=os.environ['HOME']+'/NDM/bin/rundm90'
   #eampot=os.environ['HOME']+"/Flow/INPUT/eamtab.potin_M2EX_25p_new"
   eampot=os.environ['HOME']+"/Flow/INPUT/eamtab.potin_M10_25p_new"
 if computer=='my_computer':
   ELASTIC=os.environ['HOME']+'/Codes/Elastic/'
   #EXENDM=os.environ['HOME']+'/Codes/NDM/bin/rundm90'
   EXENDM=globalv.ndm_exe
   #EXENDM=os.environ['HOME']+'/Codes/NDM/bin/rundm90'
   #eampot=os.environ['HOME']+"/Codes/Flow/INPUT/eamtab.potin_M2EX_25p_new"
 if computer=='irene':
   ELASTIC=os.environ['HOME']+'/FitFormation/Elastic/'
   EXENDM=os.environ['HOME']+'/NDM/bin/rundm90'
   #eampot=os.environ['HOME']+"/Flow/INPUT/eamtab.potin_M2EX_25p_new"
   eampot=os.environ['HOME']+"/Flow/INPUT/eamtab.potin_M10_25p_new"

 eampot=globalv.eampot
 nodes=globalv.nodes
 ntasks=globalv.ntasks
 root_dir=globalv.root_dir
 dirNDM=root_dir+'/NDM'
 dinfile=globalv.input_din_min
 mlfile=globalv.input_ml
 # relaxation = 'yes' / 'no'
 # if the relaxation is performed.E.g. no - in the case of bulk structures; yes - for SIAs or others unrelaxed strunctures

 #os.system('pwd')

 if not os.path.exists(dinfile):
    print("ndm setup_ndm: DIN file %s doesn't exist. Put the correct path in setup_ndm"%str(dinfile))
    exit(0)


 if not os.path.exists(mlfile):
    print("ndm setup_ndm: ML file %s doesn't exist. Put the correct path in setup_ndm"%str(mlfile))
    exit(0)


 if not os.path.exists(eampot):
    print("ndm setup_ndm: Potential file %s doesn't exist. Put the correct path in setup_ndm"%str(eampot))
    exit(0)

 cal_elas=globalv.cal_elas
 if (cal_elas=='qha_bulk') or (cal_elas=='qha_cxx'):
    #limit for zero in frequencies (THz untis).  In order to detect  the three PBC modes
    limit_omega=globalv.limit_omega
    relaxation=globalv.relaxation
    #should be nodes multiply by 24 ntasks-per-node
    nodes=globalv.nodes
    ntasks=globalv.ntasks
    eigenvalues_phondy="eigenvalues.dat"
    #dinfile_ph=dirNDM+'/PHONDY/cube.din_ph'
    dinfile_ph=globalv.input_din_ph
    if computer=='zelda':
      EXENDM_PH=os.environ['HOME']+ '/NDM/bin/rundm90_ml_phondy_para'
      SAUVE2GIN=os.environ['HOME']+ '/bin/sauve2gin'
    if computer=='my_computer':
      EXEMIN_PH=os.environ['HOME']+'/home/marinica/Flow/lammps.git/src/lmp_mpi'
      EXENDM_PH=globalv.ndm_exe_ph
      #EXENDM_PH=os.environ['HOME']+ '/NDM/p2.8/rundm90_phondy_para_lammps'
      SAUVE2GIN=os.environ['HOME']+ '/bin/sauve2gin'
    if computer=='occigen':
      #EXENDM_PH=os.environ['HOME']+ '/NDM/bin/rundm90_ml_phondy_para'
      EXENDM_PH=globalv.ndm_exe_ph
      SAUVE2GIN=os.environ['HOME']+ '/bin/sauve2gin'
    if computer=='irene':
      #EXENDM_PH=os.environ['HOME']+ '/NDM/bin/rundm90_ml_phondy_para'
      EXEMIN_PH=globalv.minexe_ph
      EXENDM_PH=globalv.ndm_exe_ph
      SAUVE2GIN=os.environ['HOME']+ '/bin/sauve2gin'
      
      
      
    if not os.path.exists(dinfile_ph):
       print("ndm setup_ndm_ph: DIN PH file %s doesn't exist. \
          Put the correct path in setup_ndm_ph"%str(dinfile_ph))
       exit(0)
    flmodes=dirNDM+'/PHONDY/cube.phondy.lmodes'
    fldos  =dirNDM+'/PHONDY/cube.phondy.ldos'
    inpmab =dirNDM+'/PHONDY/cube.phondy'

    #dinfile_cg=dirNDM+'/PHONDY/cube.din_cg'
    dinfile_cg=globalv.input_din_min
    dinfile_cg2=dirNDM+'/PHONDY/cube.din_cg2'
    dinfile_tr=dirNDM+'/PHONDY/cube.din_tr'


    for testf in (EXENDM_PH, dinfile_ph, inpmab, flmodes, fldos, dinfile_cg, dinfile_cg2,dinfile_tr):
      if not os.path.exists(testf):
        print("ndm setup_ndm_ph 1: The file %s doesn't exist. Put the correct path. "%(str(testf)))
        exit(0)


    minexe_ph=globalv.minexe_ph
    EXEMIN_PH=minexe_ph
    mininput_ph=globalv.mininput_ph

    minpot01_ph=globalv.minpot01_ph
    minpot02_ph=globalv.minpot02_ph
    for testf in (minexe_ph, mininput_ph, minpot01_ph, minpot02_ph):
      if not os.path.exists(testf):
        print("ndm setup_ndm_ph 2: The file %s doesn't exist. Put the correct path. "%(str(testf)))
        exit(0)


 return ;
