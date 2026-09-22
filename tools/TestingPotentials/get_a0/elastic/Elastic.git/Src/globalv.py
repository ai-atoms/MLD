import os
from lammps import run_lammps
from ndm    import run_ndm
from phondy import run_phondy
from vasp   import run_vasp, run_vasp_ph
from pwscf  import run_pwscf
#init.py

def init():
 global nt_strain
 global perc_strain_min,perc_strain_max
 global dt
 global atoms_per_unit
#input parameters defined global
 global mode, abinitio
 global typerun, debug
 global cluster, nodes, ntasks, npools,type_node
 global structure
 global prefix

 global localPWSCF, localVASP
 global input_dirVASP, input_dirLAMMPS,  input_incar, input_neb_incar, input_kpoint, input_potcar
 global dirVASP,  incar, kpoint, potcar
 global eampot, ndmexe, input_din_min, input_din_ph, input_ml
 global dynamics

# runtime  variables
 global based
 global root_dir
 global script
# elastic
 global a_guess
 global cal_elas, bulk_deform
 global number_of_deformations
 global index_list_of_deformations_nostrandard, value_list_of_deformations_nostrandard, list_of_deformations_nostandard
 global beta_matrix
 global aneto, aneto_exe
 global ndm_exe, ndm_exe_ph, lammps_exe 
#qha

 global temp_min,temp_max
 global perc_vol_min,perc_vol_max
 global nt_temp
 global limit_omega, limit_sigma, relaxation
 global nodes, ntasks, cpu_per_task, runtime, project, queue
 #active only for qha AND vasp
 global type_exe_system
 global disp_vaspi
 global nt_strain_min, nt_strain_max
 global nt_bulk_run
 global pre_build_requested
# energy, volume
 global energy,valc11,valc12,valc44
 global volume
 global nat_per_box
 global pack_disp_size

#input parameters defined global
 global a0base
#global conversion
 global unit_conversion1
 global minexe_ph, mininput_ph, phinput_ph, minpot01_ph, minpot02_ph, qha_minmode, qha_phmode
 global lammps_to_gin, gin_to_lammps
 
 minexe_ph=''
 mininput_ph=''
 phinput_ph=''
 minpot01_ph=''
 minpot02_ph=''
 qha_phmode=''
 qha_minmode=''
 lammps_to_gin=''
 gin_to_lammps=''

 unit_conversion1 = 160.217646 # 1 eV/A^3 = 1.60217646E-19 / 1E-30 Pa = 160.217646 GPa
 mode=''
 typerun=''
 debug=0
 structure=''
 based=''
 root_dir=''
 script=''
 dt=0.1
 prec_strain_min=0.0
 prec_strain_max=0.0
 atoms_per_unit=1
 energy=[]
 volume=[]
 relaxation='no'
 dynamics='min'
 bulk_deform='no'
 aneto='no'
 aneto_exe=' '
 nt_bulk_run=0
 nt_strain_min=0
 nt_strain_max=1
 cpu_per_task=1
def run_md(mode, based, dirrun,Exceptions=None):

 os.chdir(based+'/'+dirrun)
 test=0
 if Exceptions is None:
    test=1
 else:
    for tmp in Exceptions:
       if str(tmp)==str(mode):
         print(('mode %s is not implemented'%mode))
         exit(0)
       else:
         test=1


 if test==1 :
   if ( mode == 'ndm') :
     run_ndm("out.run")

   if ( mode == 'phondy'):
     run_phondy("out.run")


   if ( mode == 'vasp') :
     if (cal_elas == 'qha_bulk'):
        print(dirrun, nt_bulk_run, based, dirrun)
        if  dirrun=='bulk':
           if  (nt_bulk_run==0) :
              print('no running in bulk')
           else:
              run_vasp_ph(based, dirrun, "out.run")
        else:
           run_vasp_ph(based, dirrun, "out.run")
     else:
          run_vasp("out.run")

   if ( mode == 'pwscf') :
     run_pwscf("out.run")


   if ( mode == 'lammps') :
     if (cal_elas == 'qha_bulk'):
        if  dirrun=='bulk':
           if  (nt_bulk_run==0) :
              print('no running in bulk')
           else:
              run_lammps_ph(based, dirrun, "out.run")
        else:
           run_lammps_ph(based, dirrun, "out.run")
     else:
          run_lammps(based, dirrun, "out.run")


 return;
