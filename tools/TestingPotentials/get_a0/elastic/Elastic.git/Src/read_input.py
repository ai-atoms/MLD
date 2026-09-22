import configparser
import globalv
import string
import os
import numpy as np

#class Config_perso(ConfigParser):

#  def __init__


def get_perso(self, option, default, section):
        """
        get a value from the config but return
        a default value if it is not found
        The type of default defines the type
        returned
        """
        print('debug', section)
        try:
            if type(default) == int:
                res = self.getint(option, section)
            elif type(default) == float:
                res = self.getfloat(option, section)
            elif type(default) == bool:
                res = self.getboolean(option, section)
            else:
                res = self.getstr(option, section)
        except (NoOptionError, NoSectionError):
            return default
        return res



def read_input_file(input_ini) :
 #ConfigParser.NoOptionError()
 #config = ConfigParser.ConfigParser(allow_no_value=True)
 config = configparser.ConfigParser()
 config.read(input_ini)

 try:
   globalv.mode=str(config.get("Common","mode")).strip()
 except configparser.NoOptionError:
   print('The mode should be given in Common section: vasp, ndm, phondy, lammps, pwscf')
   exit(0)

 try:
   globalv.typerun=str(config.get("Common","typerun")).strip()
 except configparser.NoOptionError:
   print('The typerun should be given in Common section: pre_build, build, run, extract, post_extract')
   exit(0)


 globalv.cluster=str(config.get("Common","cluster")).strip()

 try:
   globalv.dynamics=str(config.get("Common","dynamics")).strip()
 except configparser.NoOptionError:
   globalv.dynamics=str('min')

 try:
   globalv.debug=int(config.get("Common","debug"))
 except configparser.NoOptionError:
   globalv.debug=0



 try:
   globalv.abinitio=int(config.get("Common","abinitio"))
 except configparser.NoOptionError:
   globalv.abinitio=0






 globalv.a_guess = float(config.get("Common","a_guess"))
 globalv.structure=str(config.get("Common","structure")).strip()

 globalv.nodes=int(config.get("Common","nodes"))
 globalv.ntasks=int(config.get("Common","ntasks"))
 globalv.npools=int(config.get("Common","npools"))

 globalv.prefix=str(config.get("Common","prefix")).strip()
 globalv.type_node=str(config.get("Common","type_node")).strip()
 globalv.cpu_per_task=str(config.get("Common","cpu_per_task")).strip()
 globalv.runtime=str(config.get("Common","runtime")).strip()
 globalv.project=str(config.get("Common","project")).strip()
 globalv.queue=str(config.get("Common","queue")).strip()

 globalv.perc_strain_min= float(config.get("Common","percent_for_strain_min"))/100.0
 globalv.perc_strain_max= float(config.get("Common","percent_for_strain_max"))/100.0
 globalv.nt_strain = int(config.get("Common","nt_strain"))


 if (globalv.mode == 'ndm' ) or \
    (globalv.mode == 'vasp' ) or \
    (globalv.mode == 'pwscf' ) or \
    (globalv.mode == 'lammps' ) or \
    (globalv.mode == 'phondy' ):
    globalv.cal_elas=str(config.get("Common","cal_elas")).strip()


 if globalv.cal_elas=='polar':
    try:
      globalv.bulk_deform=str(config.get("Common","bulk_deform")).strip()
    except configparser.NoOptionError:
      globalv.bulk_deform=str('no')

 try:
   globalv.aneto=str(config.get("Common","aneto")).strip()
   if not ((globalv.aneto=='yes') or (globalv.aneto=='no')):
      print('aneto option accept only yes ot no. Stop')
      exit(0)
 except configparser.NoOptionError:
   globalv.aneto=str('no')

 if globalv.aneto=='yes':
    globalv.aneto_exe=os.environ['HOME']+'/'+'/bin/aneto.x'
    if not os.path.exists(globalv.aneto_exe):
     print("read_input: Aneto exe %s doesn't exist. "%str(aneto_exe))
     exit(0)


 if (globalv.mode=='vasp'):

    globalv.input_dirVASP=str(config.get("Vasp","input_dirVASP")).strip()
    globalv.input_incar=str(config.get("Vasp","input_incar")).strip()
    globalv.input_kpoints=str(config.get("Vasp","input_kpoints")).strip()
    globalv.input_potcar=str(config.get("Vasp","input_potcar")).strip()


 if (globalv.mode=='lammps'):
    try:
       globalv.input_dirLAMMPS=str(config.get("LAMMPS","input_dirLAMMPS")).strip()
    except configparser.NoOptionError:
      globalv.input_dirLAMMPS=str('LAMMPS')

    globalv.input_incar=str(config.get("LAMMPS","input_incar")).strip()
    try:
      globalv.input_neb_incar=str(config.get("LAMMPS","input_neb_incar")).strip()
    except  configparser.NoOptionError:
      globalv.input_neb_incar=str('lammps_neb_fire.in')
      
      
    try:
      globalv.lammps_exe=str(config.get("LAMMPS","lammps_exe")).strip()
    except  configparser.NoOptionError:
      globalv.lammps_exe=str('lmp_mpi')

    #globalv.input_kpoints=str(config.get("LAMMPS","EXE_file")).strip()
    globalv.input_potcar=str(config.get("LAMMPS","input_potcar")).strip()





 if (globalv.mode=='ndm'):
    try:
        globalv.eampot=str(config.get("NDM","eampot")).strip()
    except  configparser.NoOptionError:
        globalv.eampot=str('/home/marinica/Flow/INPUT/eamtab.potin_M2EX_25p_new')

    try:
        globalv.input_din_min=str(config.get("NDM","input_din_min")).strip()
    except  configparser.NoOptionError:
        globalv.input_din_min=str('cube.din_min')



    try:
        globalv.input_din_ph=str(config.get("NDM","input_din_ph")).strip()
    except  configparser.NoOptionError:
        globalv.input_din_ph=str('cube.din_ph')


    try:
        globalv.input_ml=str(config.get("NDM","input_ml")).strip()
    except  configparser.NoOptionError:
        globalv.input_ml=str('cube.ml_min')


    try:
        globalv.ndm_exe=str(config.get("NDM","ndmexe")).strip()
    except  configparser.NoOptionError:
        globalv.ndm_exe=str(os.environ['HOME']+'/Codes/NDM/bin/rundm90')


    try:
        globalv.ndm_exe_ph=str(config.get("NDM","ndmexe_ph")).strip()
    except  configparser.NoOptionError:
        globalv.ndm_exe_ph=str(os.environ['HOME']+'/Codes/NDM/bin/rundm90_phondy')

 if  globalv.cal_elas=='qha_bulk' or globalv.cal_elas=='qha_cxx':
  globalv.perc_vol_min=float(config.get("Qha","percent_vol_min"))/100.0
  globalv.perc_vol_max=float(config.get("Qha","percent_vol_max"))/100.0
  globalv.temp_min=float(config.get("Qha","temp_min"))
  globalv.temp_max=float(config.get("Qha","temp_max"))
  globalv.nt_temp=int(config.get("Qha","nt_temp"))
  globalv.relaxation=str(config.get("Qha","relaxation")).strip()
  globalv.limit_omega=float(config.get("Qha","limit_omega"))


  try:
    globalv.minexe_ph=str(config.get("Qha","minexe_ph")).strip()
  except  configparser.NoOptionError:
    globalv.minexe_ph=str(os.environ['HOME']+'/Codes/NDM/bin/rundm90')

  try:
    globalv.mininput_ph=str(config.get("Qha","mininput_ph")).strip()
  except  configparser.NoOptionError:
    globalv.mininput_ph=str(os.environ['HOME']+'/Codes/NDM/bin/rundm90')


  try:
    globalv.phinput_ph=str(config.get("Qha","phinput_ph")).strip()
  except  configparser.NoOptionError:
    globalv.phinput_ph=str(os.environ['HOME']+'/Codes/NDM/bin/rundm90')


  try:
    globalv.minpot01_ph=str(config.get("Qha","minpot01_ph")).strip()
  except  configparser.NoOptionError:
    globalv.minpot01_ph=str(os.environ['HOME']+'/Codes/NDM/bin/rundm90')

  try:
    globalv.minpot02_ph=str(config.get("Qha","minpot02_ph")).strip()
  except  configparser.NoOptionError:
    globalv.minpot02_ph=str(os.environ['HOME']+'/Codes/NDM/bin/rundm90')

  try:
    globalv.qha_minmode=str(config.get("Qha","qha_minmode")).strip()
  except  configparser.NoOptionError:
    globalv.qha_minmode=str('ndm')


  try: 
    globalv.lammps_to_gin=str(config.get("Qha","lammps_to_gin")).strip()
  except: 
    globalv.lammps_to_gin=str('lammps_to_gin.py')

  try: 
    globalv.gin_to_lammps=str(config.get("Qha","gin_to_lammps")).strip()
  except: 
    globalv.lammps_to_gin=str('gin_to_lammps.py')

      
  try:
    globalv.qha_phmode=str(config.get("Qha","qha_phmode")).strip()
  except  configparser.NoOptionError:
    globalv.qha_phmode=str('phondy-ndm')

  if globalv.mode=='vasp':
    globalv.type_exe_system=str(config.get("Qha","type_exe_system")).strip()
    globalv.list_of_disp=list(map(int, config.get("Qha","list_of_disp").split()))

    try:
       globalv.nt_strain_min=int(config.get("Qha","nt_strain_min"))
    except  configparser.NoOptionError:
       globalv.nt_strain_min=int(0)
    try:
       globalv.nt_strain_max=int(config.get("Qha","nt_strain_max"))
    except  configparser.NoOptionError:
       globalv.nt_strain_max=int(globalv.nt_strain-1)


    try:
       globalv.pack_disp_size=int(config.get("Qha","pack_disp_size"))
    except  configparser.NoOptionError:
       globalv.pack_disp_size=int(1)
    if globalv.pack_disp_size > 1 and not(globalv.cluster == 'occigen' or globalv.cluster == 'marconi' or globalv.cluster=='irene') :
         print('read_input ... pack_disp_size > 1 implemented only for occigen / marconi / irene')
         exit(0)

    globalv.nt_bulk_run=int(config.get("Qha","nt_bulk_run"))


 if globalv.script=='elastic.py':
   #globalv.perc= float(config.get("Elastic","percent_for_strain"))

   try:
      globalv.limit_sigma = float(config.get("Elastic","limit_sigma"))
   except configparser.NoOptionError:
      globalv.limit_sigma = float(1.e-5)

 if globalv.cal_elas=='polar':
    #globalv.number_of_deformations= float(config.get("Elastic","number_of_deformations"))
    globalv.list_of_deformations=list(map(int,config.get("Elastic","list_of_deformations").split()))
    globalv.valc11=float(config.get("Elastic","c11"))
    globalv.valc12=float(config.get("Elastic","c12"))
    globalv.valc44=float(config.get("Elastic","c44"))


    try:
      index_list_of_deformations_nostandard=[]
      value_list_of_deformations_nostandard=[]
      globalv.list_of_deformations_nostandard=list(map(float,config.get("Elastic","list_of_deformations_nostandard").split()))
      tmp=np.array(globalv.list_of_deformations_nostandard)
      if len(tmp)%2!=0:
        print("error in input of list_of_deformations_nostandard")
        exit(0)
      else:
        for i in range(len(tmp)):
          if i%2==0:
            index_list_of_deformations_nostandard.append(int(tmp[i]))
          if i%2==1:
            value_list_of_deformations_nostandard.append(float(tmp[i])/100.0)

      globalv.index_list_of_deformations_nostandard=index_list_of_deformations_nostandard
      globalv.value_list_of_deformations_nostandard=value_list_of_deformations_nostandard

    except  configparser.NoOptionError:
      globalv.list_of_deformation_nostandard=""






    if (globalv.structure != "gin"):
      print("cal_elas moe=polar is compatible only with structure=gin")
      print("now structure  = %s"%structure)
      exit(0)



 if globalv.cal_elas=='qha_cxx':
    #globalv.number_of_deformations= float(config.get("Elastic","number_of_deformations"))
    try:
      globalv.list_of_deformations=list(map(int,config.get("Elastic","list_of_deformations").split()))
    except configparser.NoOptionError:
      globalv.list_of_deformations=list(map(int," 1 2 3 ".split()))



 return ;
