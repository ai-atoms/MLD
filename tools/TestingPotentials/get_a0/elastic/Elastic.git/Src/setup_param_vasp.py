import globalv
import os

def init():
 global incar, kpoints, potcar, dirVASP,prefix
 global VASPX
 global nodes, ntasks,nkpar,type_node
 global EXENDM_PH, dinfile_ph_01, dinfile_ph_02, PHON2PREPA
 global list_of_disp
 global flmodes, fldos,eigenvalues_phondy

 base=os.getcwd()

 # occigen, zelda, my_computer
 computer=globalv.cluster

 if computer=='occigen':
   ELASTIC=os.environ['HOME']+'/Elastic/'
   VASPX=os.environ['HOME']+'/Vasp5.4/vasp.5.4.1/bin/vasp_std'
   #VASPX=os.environ['STOREDIR']+'/wPH/NDM/bin/rundm90'
 if computer=='marconi':
   ELASTIC=os.environ['PRA']+'/Elastic/'
   VASPX='/cineca/prod/opt/applications/vasp/5.4.1_knl/intelmpi--2017--binary/bin/vasp-std'
   #VASPX=os.environ['HOME']+'/Vasp5.4/vasp.5.4.1/bin/vasp_std'
   #VASPX=os.environ['STOREDIR']+'/wPH/NDM/bin/rundm90'
 if computer=='zelda':
   ELASTIC=os.environ['HOME']+'/Elastic/'
   VASPX=os.environ['HOME']+'/Vasp5.4/vasp.5.4.1/bin/vasp_std'
 if computer=='my_computer':
   ELASTIC=os.environ['HOME']+'/Elastic/'
   VASPX=os.environ['HOME']+'/Vasp5.4/vasp.5.4.1/bin/vasp_std'
 if computer=='irene':
  try:
   ELASTIC=os.environ['CCCWORKDIR']+'/Flow/Elastic.git/'
  except KeyError:
   ELASTIC=os.environ['HOME']+'/Elastic.git/'     

  VASPX='/ccc/products/vasp-6.4.0/intel--20.0.0__openmpi--4.0.1/default/bin/vasp_std'

 root_dir=globalv.root_dir
 dirVASP=root_dir+"/"+globalv.input_dirVASP
 incar=dirVASP+"/"+globalv.input_incar
 kpoints=dirVASP+"/"+globalv.input_kpoints
 potcar=dirVASP+"/"+globalv.input_potcar
 prefix=globalv.prefix
 nodes=globalv.nodes
 ntasks=globalv.ntasks 
 nkpar=globalv.npools
 type_node=globalv.type_node

 # relaxation = 'yes' / 'no'
 # if the relaxation is performed.E.g. no - in the case of bulk structures; yes - for SIAs or others unrelaxed strunctures
 

 if not os.path.exists(incar):
    print("vasp setup_vasp: incar  file %s doesn't exist. Put the correct path"%str(incar))
    exit(0)

 if not os.path.exists(kpoints):
    print("vasp setup_vasp: kpoints  file %s doesn't exist. Put the correct path"%str(kpoints))
    exit(0)
 if not os.path.exists(potcar):
    print("vasp setup_vasp:potcar   file %s doesn't exist. Put the correct path"%str(potcar))
    exit(0)


 script=globalv.script
 cal_elas=globalv.cal_elas
 if (cal_elas=='qha_bulk' or cal_elas=='qha_cxx'):
    
    #limit for zero in frequencies (THz untis).  In order to detect  the three PBC modes
    limit_omega=globalv.limit_omega
    relaxation=globalv.relaxation
    list_of_disp=globalv.list_of_disp
    #should be nodes multiply by 24 ntasks-per-node
    eigenvalues_phondy="eigenvalues.dat"

    dinfile_ph_01=dirVASP+'/PHONDY/INPHON_01'
    dinfile_ph_02=dirVASP+'/PHONDY/INPHON_02'
    if computer=='zelda':
      EXENDM_PH=os.environ['HOME']+ '/NDM/p2.8/rundm90_phondy_para'
      EXENDM_PH_B=os.environ['HOME']+ '/Elastic/Zelda/Bin/phon_bulk.x'
      EXENDM_PH_D=os.environ['HOME']+ '/Elastic/Zelda/Bin/phon_defect.x'
      EXENDM=os.environ['HOME']+ '/NDM/p2.8/rundm90_phondy_para'
      SAUVE2GIN=os.environ['HOME']+ '/bin/sauve2gin'
    if computer=='my_computer':
      EXENDM_PH_B=os.environ['HOME']+ '/Elastic/Bin/MyComputer/phon_bulk.x'
      EXENDM_PH_D=os.environ['HOME']+ '/Elastic/Bin/MyComputer/phon_defect.x'
      PHON2PREPA=os.environ['HOME']+ '/Elastic/Bin/MyComputer/phon2prepa.x'
      SAUVE2GIN=os.environ['HOME']+ '/bin/sauve2gin'

    if computer=='occigen':
      EXENDM_PH_B=os.environ['HOME']+ '/Elastic/Bin/Occigen/phon_bulk.x'
      EXENDM_PH_D=os.environ['HOME']+ '/Elastic/Bin/Occigen/phon_defect.x'
      PHON2PREPA=os.environ['HOME']+ '/Elastic/Bin/Occigen/phon2prepa.x'
      SAUVE2GIN=os.environ['HOME']+ '/bin/sauve2gin'
    if computer=='marconi':
      EXENDM_PH_B=os.environ['PRA']+ '/Elastic/Bin/Occigen/phon_bulk.x'
      EXENDM_PH_D=os.environ['PRA']+ '/Elastic/Bin/Occigen/phon_defect.x'
      PHON2PREPA=os.environ['PRA']+ '/Elastic/Bin/Occigen/phon2prepa.x'
      SAUVE2GIN=os.environ['PRA']+ '/bin/sauve2gin'
    if computer=='irene':
      EXENDM_PH_B=ELASTIC + '/Bin/Occigen/phon_bulk.x'
      EXENDM_PH_D=ELASTIC + '/Bin/Occigen/phon_defect.x'
      PHON2PREPA=ELASTIC + '/Bin/Occigen/phon2prepa.x'
      SAUVE2GIN=os.environ['HOME']+ '/bin/sauve2gin'


      #EXENDM_PH=os.environ['STOREDIR']+ '/wPH/NDM/bin/rundm90_phondy_para'
      #SAUVE2GIN=os.environ['HOME']+ '/bin/sauve2gin'

    if globalv.type_exe_system=='defect':
        EXENDM_PH=EXENDM_PH_D
    elif globalv.type_exe_system=='bulk':
        EXENDM_PH=EXENDM_PH_B
    else:
         print("the only choice for type_exe_system is bulk or defects")
         print(("current value %s"%type_exe_system))
         exit(0)
    # this will probably should be removed in the future
    flmodes=dirVASP+'/PHONDY/cube.phondy.lmodes'
    fldos  =dirVASP+'/PHONDY/cube.phondy.ldos'
    inpmab =dirVASP+'/PHONDY/cube.phondy'

    dinfile_cg=dirVASP+'/PHONDY/cube.din_cg'
    dinfile_cg2=dirVASP+'/PHONDY/cube.din_cg2'
    dinfile_tr=dirVASP+'/PHONDY/cube.din_tr'


    #for testf in (EXENDM_PH, dinfile_ph, inpmab, flmodes, fldos, dinfile_cg, dinfile_cg2,dinfile_tr):
    for testf in (EXENDM_PH, dinfile_ph_01, dinfile_ph_02, PHON2PREPA)  :
      if not os.path.exists(testf):
        print("vasp setup_param_vasp: The file %s doesn't exist. Put the correct path. "%(str(testf)))
        exit(0)

 return ;
