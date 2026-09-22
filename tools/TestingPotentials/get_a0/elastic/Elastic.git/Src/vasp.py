import globalv
import os
import re
import sys
import math
import setup_param_vasp
import numpy as np
import mmap
from tools import search_word_replace_line


def setup_vasp(based,indir,cname):

# I put a commnet
  os.chdir(based+'/'+indir)
  setup_param_vasp.init()
  incar=setup_param_vasp.incar
  kpoints=setup_param_vasp.kpoints
  potcar=setup_param_vasp.potcar
  nkpar=setup_param_vasp.nkpar
  computer=globalv.cluster
  typerun=globalv.typerun
  prefix=setup_param_vasp.prefix

  os.system("cp -f %s  INCAR"%(incar))
  os.system("ln -s %s  POTCAR"%(potcar))
  os.system("ln -s %s  KPOINTS"%(kpoints))

  if typerun=='pre_build':
     if globalv.relaxation=='yes':
       search_word_replace_line('INCAR', 'NSW', ' NSW =200' +'\n')

     if globalv.relaxation=='no':
       search_word_replace_line('INCAR', 'NSW', ' NSW =1' +'\n')


  if computer=='zelda':
     generate_jsub_vasp_in_zelda(indir,cname)
  if computer=='occigen':
     generate_jsub_vasp_in_occigen(indir,cname)
  if computer=='irene':
     generate_jsub_vasp_in_irene(indir,cname)
  if computer=='marconi':
     generate_jsub_vasp_in_marconi(indir,cname)




  return;


def setup_vasp_ph(based,indir,cname):

  os.chdir(based+'/'+indir)
  setup_param_vasp.init()
  incar=setup_param_vasp.incar
  kpoints=setup_param_vasp.kpoints
  potcar=setup_param_vasp.potcar
  nkpar=setup_param_vasp.nkpar
  computer=globalv.cluster
  prefix=setup_param_vasp.prefix
  list_of_disp=setup_param_vasp.list_of_disp
  dinfile_ph_01=setup_param_vasp.dinfile_ph_01
  dinfile_ph_02=setup_param_vasp.dinfile_ph_02
  EXENDM_PH=setup_param_vasp.EXENDM_PH
  PHON2PREPA=setup_param_vasp.PHON2PREPA
  #if not os.path.exists('min/CONTCAR'):
  #      print "Probably you need a relaxation "%(str(testf))
  #      exit(0)


  #os.system("cp  min/CONTCAR POSCAR")

  for i in range(len(list_of_disp)):
     os.chdir(based+'/'+indir)
     disp=list_of_disp[i]
     dir_mresults='mresults_'+str(int(disp))
     os.mkdir(dir_mresults)
     os.chdir(dir_mresults)
     os.system('cp %s  INPHON'%(dinfile_ph_01))
     if globalv.relaxation=='no':
      os.system('cp %s/POSCAR   POSCAR'%(based+'/pre'+indir))
      os.system('cp %s/POSCAR   OLDPOSCAR'%(based+'/pre'+indir))
     if globalv.relaxation=='yes':
      os.system('cp %s/CONTCAR   CONTCAR'%(based+'/pre'+indir))
      if not os.path.exists(based+'/pre'+indir+'/CONTCAR'):
        print("CONTCAR absent in %s"%(based+'/pre'+indir))
      else :
        with open ("CONTCAR","r") as textobj:
          listf=list(textobj)
        no_line=6
        del listf[no_line-1]
        with open("POSCAR","w") as textnew:
          for n in listf:
             textnew.write(n)
        os.system('cp %s/CONTCAR   OLDCONTCAR'%(based+'/pre'+indir))

     search_word_replace_line('INPHON', 'DISP', 'DISP='+str(disp) +'\n')
     os.system('%s'%EXENDM_PH)
     os.system('mkdir -p tmpstr')
     os.system('%s POSCAR l.gin'%PHON2PREPA)
     os.system('mv *.poscar tmpstr')
     #os.system('mv *.pwscf tmpstr')
     os.system('rm -f  *.pwscf ')
     os.system('rm -f  *.phondy')
     os.system('mv *.inp tmpstr')
     os.system('mv list  tmpstr')
     with open('tmpstr/list') as flist:
       fdisp=flist.read().splitlines()
     for tmp in fdisp:
       tmp_dir=tmp.strip()
       os.system('rm -rf %s'%tmp_dir)
       os.mkdir(tmp_dir)
       os.system('head -7 POSCAR > %s/POSCAR'%tmp_dir)
       os.system('cat tmpstr/%s.poscar >> %s/POSCAR'%(tmp_dir,tmp_dir))
       os.system('cp %s %s/INCAR'%(incar,tmp_dir))
       os.chdir(tmp_dir)


       if computer=='zelda':
         generate_jsub_vasp_in_zelda('M'+indir+'m'+str(disp),'d'+str(tmp_dir))
       if computer=='occigen':
         generate_jsub_vasp_in_occigen('M'+indir+'m'+str(disp),'d'+str(tmp_dir))
       if computer=='marconi':
         generate_jsub_vasp_in_marconi('M'+indir+'m'+str(disp),'d'+str(tmp_dir))
       if computer=='irene':
         generate_jsub_vasp_in_irene('M'+indir+'m'+str(disp),'d'+str(tmp_dir))



       os.system("ln -s %s  POTCAR"%(potcar))
       os.system("ln -s %s  KPOINTS"%(kpoints))

       #search_word_replace_line('INCAR', 'IBRION', 'IBRION=-1' +'\n')
       os.chdir('../')
     if globalv.pack_disp_size > 1:
      icount=0
      vlabel=100000
      ilabel=0
      for tmp in fdisp:
        tmp_dir=tmp.strip()
        if icount%globalv.pack_disp_size==0:
          if computer=='occigen':
            ilabel=ilabel+1
            clabel=str(vlabel+ilabel)
            generate_jsub_vasp_in_occigen('M'+indir+'m'+str(disp),'d'+str(tmp_dir),IdFile=clabel[-3:],NoLastPart=True)
          if computer=='marconi':
            ilabel=ilabel+1
            clabel=str(vlabel+ilabel)
            generate_jsub_vasp_in_marconi('M'+indir+'m'+str(disp),'d'+str(tmp_dir),IdFile=clabel[-3:],NoLastPart=True)
          if computer=='irene':
            ilabel=ilabel+1
            clabel=str(vlabel+ilabel)
            generate_jsub_vasp_in_irene('M'+indir+'m'+str(disp),'d'+str(tmp_dir),IdFile=clabel[-3:],NoLastPart=True)
            
        f=open('jsub'+'_'+clabel[-3:],'a')
        f.write("cp -rp ${SUBMISSION}/%s ${VASP_TMPDIR}/%s\n"%(tmp_dir, tmp_dir))
        f.write("cd     ${VASP_TMPDIR}/%s\n"%(tmp_dir))
        if (computer=="irene"):
          f.write("ccc_mprun  ${EXE}   |tee ${SUBMISSION}/wfile.out \n")
        else:   
          f.write("srun  --mpi=pmi2 -K1 --resv-ports -n ${SLURM_NTASKS}  ${EXE}   |tee ${SUBMISSION}/wfile.out \n")
        f.write("cp OUTCAR  ${SUBMISSION}/%s/\n"%(tmp_dir))
        f.write("cp CONTCAR ${SUBMISSION}/%s/\n"%(tmp_dir))
        f.write("cp OSZICAR ${SUBMISSION}/%s/\n"%tmp_dir)
        #f.write("cd .. \n")
        f.write("\n")
        icount=icount+1

  return;






def run_vasp(outVASP):

 script=globalv.script
 setup_param_vasp.init()
 computer=globalv.cluster


 os.system('pwd')
# fj=open('jsub','a')
# fj.write('mv fe.scf.out %s'%outVASP)
# fj.close()


 if computer=='zelda':
   run_vasp_zelda(outVASP)

 if computer=='occigen':
   run_vasp_occigen(outVASP)

 if computer=='marconi':
   run_vasp_marconi(outVASP)

 if computer=='irene':
   run_vasp_irene(outVASP)

 if computer=='my_computer':
   print("Not yet implemented in my_computer mode")
   run_vasp_my_computer(outVASP)
 return



def run_vasp_ph(based, indir, outVASP):

 script=globalv.script
 setup_param_vasp.init()
 computer=globalv.cluster
 list_of_disp=setup_param_vasp.list_of_disp


 for i in range(len(list_of_disp)):
    os.chdir(based+'/'+indir)
    os.system('pwd')
    disp=list_of_disp[i]
    dir_mresults='mresults_'+str(int(disp))
    os.chdir(dir_mresults)

    with open('tmpstr/list') as flist:
      fdisp=flist.read().splitlines()
      print('fdisp is here: %s on %s'%(fdisp, computer))

    if globalv.pack_disp_size == 1:
       for tmp in fdisp:
         tmp_dir=tmp.strip()
         os.chdir(tmp_dir)

         if computer=='zelda':
            run_vasp_zelda(outVASP)

         if computer=='occigen':
            run_vasp_occigen(outVASP)

         if computer=='marconi':
            run_vasp_marconi(outVASP)

         if computer=='irene':
            run_vasp_irene(outVASP)

         if computer=='my_computer':
            run_ndm_my_computer(outVASP)

         os.chdir('../')

    if globalv.pack_disp_size > 1:
       icount=0
       vlabel=100000
       ilabel=0
       for tmp in fdisp:
         tmp_dir=tmp.strip()
         if icount%globalv.pack_disp_size==0:
            ilabel=ilabel+1
            clabel=str(vlabel+ilabel)
            if computer=='occigen':
               run_vasp_occigen(outVASP,IdFile=clabel[-3:])
            if computer=='marconi':
               run_vasp_marconi(outVASP,IdFile=clabel[-3:])
            if computer=='irene':
               run_vasp_irene(outVASP,IdFile=clabel[-3:])


         icount=icount+1





 return



def run_vasp_zelda(outNDM):

#  if not os.path.exists(exerun):
#   print "Exe file %s doesn't exist. Put the correct path in setup_vasp"%str(exerun)

  os.system('pwd')
  os.system('qsub jsub')
  return

def run_vasp_occigen(outNDM, IdFile=None):

#  if not os.path.exists(exerun):
#   print "Exe file %s doesn't exist. Put the correct path in setup_vasp"%str(exerun)

  os.system('pwd')
  if IdFile==None:
     os.system('sbatch  jsub')
  else:
     os.system('sbatch jsub_%s'%IdFile)
  return

def run_vasp_marconi(outNDM, IdFile=None):

#  if not os.path.exists(exerun):
#   print "Exe file %s doesn't exist. Put the correct path in setup_vasp"%str(exerun)

  os.system('pwd')
  if IdFile==None:
     os.system('sbatch  jsub')
  else:
     os.system('sbatch jsub_%s'%IdFile)
  return




def run_vasp_irene(outNDM, IdFile=None):

#  if not os.path.exists(exerun):
#   print "Exe file %s doesn't exist. Put the correct path in setup_vasp"%str(exerun)

  os.system('pwd')
  if (IdFile==None):
     os.system('ccc_msub jsub')
  else:    
     os.system('ccc_msub jsub_%s'%IdFile)
  return


def generate_jsub_vasp_in_occigen(indir,cname,IdFile=None, NoLastPart=None):

   setup_param_vasp.init()
   script=globalv.script
   cal_elas=globalv.cal_elas
   VASPX=setup_param_vasp.VASPX
   potcar=setup_param_vasp.potcar
   incar=setup_param_vasp.incar
   kpoints=setup_param_vasp.kpoints
   nodes=setup_param_vasp.nodes
   nprocs=setup_param_vasp.ntasks
   nkpar=setup_param_vasp.nkpar
   type_node=setup_param_vasp.type_node
   prefix=globalv.prefix
   if (type_node=='HSW24'):
     task_per_node=24
   elif (type_node=='BDW28'):
     task_per_node=28
   else:
     print("Unknown type of node. For occigen it should be used only:")
     print("        HSW24 or BRW28")
     exit(0)

   if (nprocs < 0) and (nodes > 0):
     nprocs=int(nodes)*task_per_node

   if (nodes < 0) and (nprocs > 0) :
     nodes = int(float(nprocs) / float(task_per_node)) + 1

   if (nodes < 0) and (nprocs < 0):
     print("At least one from nodes and nprocs should be greater than zero")
     print("nodes:  %s"%nodes)
     print("nprocs:  %s"%nprocs)
     exit(0)
   if IdFile==None:
     f=open('jsub','w')
   else:
     f=open('jsub'+'_'+IdFile,'w')
   #debug print indir, cname
   # cname - defomation003, indir is the last directory
   f.write('#!/bin/sh\n')
   f.write('#SBATCH -J c%sd%s\n'%(indir[-2:],cname[-2:]))
   f.write('#SBATCH --nodes=%s                # debug:  number of nodes \n'%nodes)
   f.write('#SBATCH --constraint=%s          # debug:  type of nodes \n'%type_node)
   f.write('#SBATCH --ntasks-per-node=%s      # debug:  tasks /  nodes \n'%task_per_node)
   f.write('#SBATCH --ntasks=%s               # the number of MPI Tasks \n'%nprocs)
   f.write('#SBATCH --threads-per-core=1      #  \n')
   f.write('#SBATCH --time  23:20:55         # Time limit: HH:MM:SS format\n')
   f.write('##SBATCH -p highmem              #  By default 2000 MB are reserved. More uncomment\n')
   f.write('#SBATCH --output  ofile.out   # stdout pf the job\n')
   f.write('#SBATCH --error   efile.err   # stderr of the job\n')
   f.write('\n')
   f.write('\n')
   f.write('module purge\n')
   f.write('module load intelmpi/2017.0.098\n')
   f.write('module load intel\n')
   f.write('module load python\n')
   f.write('##---------------------\n')
   f.write('pwd\n')
   f.write('## SLURM_SUBMIT_DIR = the directory from which sbatch was invoked \n')
   f.write('## SLURM_JOB_ID  = job ID of the executing job \n')
   f.write('\n')
   f.write('module list\n')
   f.write('\n')
   f.write('suff=%s\n'%(cname))
   f.write('export nprocs=${SLURM_NTASKS}\n')
   f.write('echo ${nprocs} \n')
   f.write('\n')
   f.write('#export OMP_NUM_THREADS=1\n')
   f.write('export EXE=%s\n'%VASPX)
   f.write('export VASP_PSEUDO=%s\n'%potcar)
   f.write('export VASP_TMPDIR=${SCRATCHDIR}/TestPP/%s/${suff}/%s\n'%(prefix, indir))
   f.write('\n')
   if globalv.abinitio==0:
     f.write('rm -rf   ${VASP_TMPDIR}\n')
     f.write('mkdir -p ${VASP_TMPDIR}\n')
   f.write('cd ${VASP_TMPDIR}\n')
   f.write('\n')
   f.write('pwd\n')
   f.write('\n')
   f.write('\n')
   f.write('\n')
   f.write('SUBMISSION=${SLURM_SUBMIT_DIR}\n')
   #f.write('cp %s                   INCAR\n'%incar)
   if NoLastPart==None:
     f.write('cp ${SUBMISSION}/INCAR  INCAR\n')
     f.write('cp $VASP_PSEUDO         POTCAR\n')
     f.write('cp ${SUBMISSION}/POSCAR  .\n')
     if globalv.abinitio!=0:
       f.write('cp ${SUBMISSION}/POSCAR   ${SUBMISSION}/POSCAR_INI\n')
       f.write('cp ${SUBMISSION}/CONTCAR  POSCAR\n')
     f.write('cp %s                  KPOINTS\n'%kpoints)
     f.write('\n')
     f.write('outputfile=${SUBMISSION}/w${suff}.out\n')
     f.write('ulimit -c unlimited\n')
     f.write('#export MPI_GROUP_MAX=1024\n')
     f.write('\n')
     f.write('srun  --mpi=pmi2 -K1 --resv-ports -n ${SLURM_NTASKS}  ${EXE}   |tee $outputfile \n')
     f.write('\n')
     f.write('cp OUTCAR  ${SUBMISSION}/\n')
     f.write('cp CONTCAR ${SUBMISSION}/\n')

   if (cal_elas != 'qha_bulk'):
     f.write('cp OSZICAR ${SUBMISSION}/\n')
   return;


def generate_jsub_vasp_in_marconi(indir,cname,IdFile=None, NoLastPart=None):

   setup_param_vasp.init()
   script=globalv.script
   cal_elas=globalv.cal_elas
   VASPX=setup_param_vasp.VASPX
   potcar=setup_param_vasp.potcar
   incar=setup_param_vasp.incar
   kpoints=setup_param_vasp.kpoints
   nodes=setup_param_vasp.nodes
   nprocs=setup_param_vasp.ntasks
   nkpar=setup_param_vasp.nkpar
   type_node=setup_param_vasp.type_node
   prefix=globalv.prefix
   if (type_node=='KNL'):
     task_per_node=40
   else:
     print("Unknown type of node. For knl n it should be used only:")
     print("        KNL")
     exit(0)

   if (nprocs < 0) and (nodes > 0):
     nprocs=int(nodes)*task_per_node

   if (nodes < 0) and (nprocs > 0) :
     nodes = int(float(nprocs) / float(task_per_node)) + 1

   if (nodes < 0) and (nprocs < 0):
     print("At least one from nodes and nprocs should be greater than zero")
     print("nodes:  %s"%nodes)
     print("nprocs:  %s"%nprocs)
     exit(0)
   if IdFile==None:
     f=open('jsub','w')
   else:
     f=open('jsub'+'_'+IdFile,'w')
   #debug print indir, cname
   # cname - defomation003, indir is the last directory
   f.write('#!/bin/sh\n')
   f.write('#SBATCH -A  Pra14_3636\n')
   f.write('#SBATCH --job-name=c%sd%s\n'%(indir[-2:],cname[-2:]))
   f.write('#SBATCH --nodes=%s               # debug:  number of nodes \n'%nodes)
   f.write('#SBATCH --ntasks-per-node=%s     # debug:  tasks /  nodes \n'%task_per_node)
   f.write('#SBATCH --cpus-per-task=1        # 1 for pure MPI \n')
   f.write('##SBATCH --constraint=%s          # debug:  type of nodes \n'%type_node)
   f.write('#SBATCH --ntasks=%s               # the number of MPI Tasks \n'%nprocs)
   f.write('##SBATCH --threads-per-core=1      #  \n')
   f.write('#SBATCH --time  23:20:55         # Time limit: HH:MM:SS format\n')
   f.write('#SBATCH --partition=knl_usr_prod \n')
   f.write('#SBATCH --output  ofile.out   # stdout pf the job\n')
   f.write('#SBATCH --error   efile.err   # stderr of the job\n')
   f.write('\n')
   f.write('\n')
   f.write('module purge\n')
   f.write('module load profile/phys\n')
   f.write('module load profile/knl\n')
   f.write('module load autoload vasp/5.4.1_knl\n')
   f.write('##---------------------\n')
   f.write('pwd\n')
   f.write('## SLURM_SUBMIT_DIR = the directory from which sbatch was invoked \n')
   f.write('## SLURM_JOB_ID  = job ID of the executing job \n')
   f.write('\n')
   f.write('module list\n')
   f.write('\n')
   f.write('suff=%s\n'%(cname))
   f.write('export nprocs=${SLURM_NTASKS}\n')
   f.write('echo ${nprocs} \n')
   f.write('\n')
   f.write('#export OMP_NUM_THREADS=1\n')
   f.write('export EXE=%s\n'%VASPX)
   f.write('export VASP_PSEUDO=%s\n'%potcar)
   f.write('export VASP_TMPDIR=${CINECA_SCRATCH}/TestPP/%s/${suff}/%s\n'%(prefix, indir))
   f.write('\n')
   if globalv.abinitio==0:
     f.write('rm -rf   ${VASP_TMPDIR}\n')
     f.write('mkdir -p ${VASP_TMPDIR}\n')
   f.write('cd ${VASP_TMPDIR}\n')
   f.write('\n')
   f.write('pwd\n')
   f.write('\n')
   f.write('\n')
   f.write('\n')
   f.write('SUBMISSION=${SLURM_SUBMIT_DIR}\n')
   #f.write('cp %s                   INCAR\n'%incar)
   if NoLastPart==None:
     f.write('cp ${SUBMISSION}/INCAR  INCAR\n')
     f.write('cp $VASP_PSEUDO         POTCAR\n')
     f.write('cp ${SUBMISSION}/POSCAR  .\n')
     if globalv.abinitio!=0:
       f.write('cp ${SUBMISSION}/POSCAR   ${SUBMISSION}/POSCAR_INI\n')
       f.write('cp ${SUBMISSION}/CONTCAR  POSCAR\n')
     f.write('cp %s                  KPOINTS\n'%kpoints)
     f.write('\n')
     f.write('outputfile=${SUBMISSION}/w${suff}.out\n')
     f.write('ulimit -c unlimited\n')
     f.write('#export MPI_GROUP_MAX=1024\n')
     f.write('\n')
     f.write('srun  --mpi=pmi2 -K1 --resv-ports -n ${SLURM_NTASKS}  ${EXE}   |tee $outputfile \n')
     f.write('\n')
     f.write('cp OUTCAR  ${SUBMISSION}/\n')
     f.write('cp CONTCAR ${SUBMISSION}/\n')

   if (cal_elas != 'qha_bulk'):
     f.write('cp OSZICAR ${SUBMISSION}/\n')
   return;


def generate_jsub_vasp_in_irene(indir,cname, IdFile=None, NoLastPart=None):

  setup_param_vasp.init()
  script=globalv.script
  cal_elas=globalv.cal_elas
  VASPX=setup_param_vasp.VASPX
  potcar=setup_param_vasp.potcar
  incar=setup_param_vasp.incar
  kpoints=setup_param_vasp.kpoints
  nodes=setup_param_vasp.nodes
  nprocs=setup_param_vasp.ntasks
  type_node=setup_param_vasp.type_node
  prefix=globalv.prefix

  if (globalv.type_node=='milan'):
     task_per_node=24
  elif (globalv.type_node=='rome'):
     task_per_node=24
  elif (globalv.type_node=='skylake'):
     task_per_node=24
  else:
     print("Unknown type of node. For irene it should be used only:")
     print("        milan or rome or skylake")
     print(globalv.type_node)
     exit(0)



  if (nprocs < 0) and (nodes > 0):
    nprocs=int(nodes)*task_per_node

  if (nodes < 0) and (nprocs > 0) :
    nodes = int(float(nprocs) / float(task_per_node)) + 1

  if (nodes < 0) and (nprocs < 0):
    print("At least one from nodes and nprocs should be greater than zero")
    print("nodes:  %s"%nodes)
    print("nprocs:  %s"%nprocs)
    exit(0)
  if IdFile==None:
    f=open('jsub','w')
  else:
    f=open('jsub'+'_'+IdFile,'w')
  f.write('#!/bin/bash\n')
  
  if IdFile==None:
    f.write('#MSUB -r %s%s\n'%(cname[:3], indir[-4:]))
  else: 
    f.write('#MSUB -r %s%s\n'%(prefix, IdFile))  
  f.write('#MSUB -n  %s\n'%nprocs)
  f.write('#MSUB -c %s\n'%globalv.cpu_per_task)

  f.write('#MSUB -q %s\n'%globalv.type_node)
  if (globalv.queue=='test'):
   f.write('#MSUB -Q test\n') 
  f.write('#MSUB -T %s\n'%globalv.runtime)
  
  if IdFile==None:
    f.write('#MSUB -o ofile.out\n')
    f.write('#MSUB -e ofile.err\n')
  else:
    f.write('#MSUB -o ofile%s.out\n'%(IdFile))
    f.write('#MSUB -o ofile%s.err\n'%(IdFile))
     
  f.write('#MSUB -A %s\n'%globalv.project)
  #f.write('#MSUB -A den\n')
  f.write('#MSUB -m scratch,work,store\n')  
  
  f.write('\n')
  f.write('suff=%s\n'%(cname))
  f.write('##---------------------------\n')
  f.write('module purge\n')
  f.write('module load mpi/openmpi/4.0.5 scalapack/mkl/21.3.0\n')
  f.write('module load fortran/inteloneapi/21.4.0\n')
  f.write('module load gnu/12.2.0\n')
  f.write('module load cmake/3.22.2\n')
  f.write('module load vasp/6.4.0\n')
  f.write('##---------------------------\n')
  f.write('#unset OMP_DISPLAY_ENV\n')
  f.write('#unset OMP_NUM_THREADS\n')
  f.write('#unset KMP_AFFINITY\n')
  f.write('#unset OMP_PLACES\n')
  f.write('#unset OMP_PROC_BIN\n')
  f.write('##------------LOLO---------------\n')
  
  
  f.write('\n')
  f.write('export nprocs=${BRIDGE_MSUB_NPROC}\n')
  f.write('export OMP_NUM_THREADS=1\n')
  #f.write('#ort EXE=${WORKDIR}/VaspC/bin/vasp_t.x\n')
  f.write('export EXE=%s\n'%VASPX)
  f.write('export VASP_PSEUDO=%s\n'%potcar)
  f.write('export VASP_TMPDIR=${CCCSCRATCHDIR}/vaspTmp/%s/${suff}/%s\n'%(prefix, indir))
  f.write('\n')
  if globalv.abinitio==0:
    f.write('rm -rf   ${VASP_TMPDIR}\n')
    f.write('mkdir -p ${VASP_TMPDIR}\n')
  f.write('cd ${VASP_TMPDIR}\n')
  f.write('\n')
  f.write('pwd\n')
  f.write('\n')
  f.write('\n')
  f.write('\n')
  f.write('SUBMISSION=${BRIDGE_MSUB_PWD}\n')
  if NoLastPart==None:
    f.write('cp ${SUBMISSION}/INCAR   INCAR\n')
    f.write('cp ${SUBMISSION}/POSCAR   POSCAR\n')
    f.write('cp ${VASP_PSEUDO}         POTCAR\n')
    if globalv.abinitio!=0:
      f.write('cp ${SUBMISSION}/POSCAR   ${SUBMISSION}/POSCAR_INI\n')
      f.write('cp ${SUBMISSION}/CONTCAR  POSCAR\n')
    f.write('cp %s                  KPOINTS\n'%kpoints)
    f.write('\n')
    f.write('pwd\n')
    if IdFile==None:
     f.write('outputfile=${BRIDGE_MSUB_PWD}/w${suff}.out\n')
    else: 
     f.write('outputfile=${BRIDGE_MSUB_PWD}/w${suff}%s.out\n'%(IdFile))   
    f.write('ulimit -s unlimited\n')
    #f.write('mpirun  -np $nprocs ${EXE}   |tee $outputfile\n')
    f.write('ccc_mprun  ${EXE}   |tee $outputfile\n')
    f.write('\n')
    f.write('cp OUTCAR  ${BRIDGE_MSUB_PWD}/OUTCAR\n')
    f.write('cp CONTCAR ${BRIDGE_MSUB_PWD}/CONTCAR\n')
  if (cal_elas != 'qha_bulk'):
     f.write('cp OSZICAR ${BRIDGE_MSUB_PWD}/\n')



  f.close()
  return;


def generate_jsub_vasp_in_zelda(indir,cname):

  cal_elas=globalv.cal_elas
  nprocs=setup_param_vasp.ntasks

  f=open('jsub','w')
  f.write('#!/bin/sh\n')
  f.write('#PBS -q prod_para\n')
  f.write('#PBS -l nodes=1:ppn=%s\n'%nprocs)
  f.write('#PBS -l walltime=10:00:20\n')
  f.write('#PBS -N c%sd%s\n'%(indir[-2:],cname[-2]))
  f.write('#PBS -e error\n')
  f.write('#PBS -o output\n')
  f.write('\n')
  f.write('suff=%s%s\n'%(indir,cname))
  f.write('\n')
  #f.write('module unload openmpi/2.0.0.\n')
  #f.write('module load  mpt/2.12\n')
  f.write('\n')
  f.write('export nprocs=12\n')
  f.write('export OMP_NUM_THREADS=1\n')
  #f.write('export EXE=/home/marinica/Flow/Vasp/bin/vasp-k-points.mpt.x\n')
  f.write('export EXE=/home/marinica/Flow/Vasp/bin/vasp-k-points_openmpi.x\n')
  f.write('\n')
  f.write('ldd $EXE\n')
  f.write('SCRATCHDIR=/scratch/marinica/CalcPP/\n')
  f.write('export VASP_TMPDIR=${SCRATCHDIR}/${suff}\n')
  f.write('mkdir -p ${VASP_TMPDIR}\n')
  f.write('pwd\n')
  f.write('cd $PBS_O_WORKDIR/\n')
  f.write('cp -p POSCAR  ${VASP_TMPDIR}\n')
  f.write('cp -p POTCAR  ${VASP_TMPDIR}\n')
  f.write('cp -p INCAR   ${VASP_TMPDIR}\n')
  f.write('cp -p KPOINTS ${VASP_TMPDIR}\n')
  f.write('cd ${VASP_TMPDIR}\n')
  f.write('pwd\n')
  f.write('ls\n')
  f.write('inputfile=$PBS_O_WORKDIR/${suff}.in\n')
  f.write('outputfile=$PBS_O_WORKDIR/w${suff}.out\n')
  f.write('ulimit -s unlimited\n')
  f.write('mpirun -n $nprocs ${EXE}   > ${outputfile}\n')
  f.write('cp -rp ${VASP_TMPDIR}/OUTCAR  $PBS_O_WORKDIR/\n')
  f.write('cp -rp ${VASP_TMPDIR}/CONTCAR  $PBS_O_WORKDIR/\n')
  f.write('rm -rf ${VASP_TMPDIR}\n')
  if (cal_elas != 'qha_bulk'):
     f.write('cp -rp ${VASP_TMPDIR}/OSZINCAR  $PBS_O_WORKDIR/\n')

  f.close()



  return


def get_natom_from_vasp_outcar_old(outVASP):
 f=open(outVASP,'r')
 for line in f:
   if re.match("(.*)number of ions(.*)", line):
     word=line.split()
     nat=float(word[11])
     break
 f.close()
 return nat;

def get_natom_from_vasp_outcar(outVASP):
 with open(outVASP,'r') as f:
   m=mmap.mmap(f.fileno(),0,prot=mmap.PROT_READ)
   i=m.find(b'number of ions')
   m.seek(i)
   line=m.readline()
   word=line.split()
   nat=float(word[-1])
 f.close()
 return nat;






def get_energy_from_vasp_outcar_old(outVASP):
 f=open(outVASP,'r')
 for line in f:
   if re.match("(.*)sigma(.*)", line):
     word=line.split()
     #print line
 f.close()
 return float(word[6]);


def get_energy_from_vasp_outcar(outVASP):
 with open(outVASP,'r') as f:
  m=mmap.mmap(f.fileno(),0,prot=mmap.PROT_READ)
  i = m.rfind(b'sigma')
  m.seek(i)
  line=m.readline()
  word=line.split()
 f.close()
 return float(word[-1]);






def get_sigma_from_vasp_old(outVASP, Patched=None):
 """
 Get sigma from patched output of VASP
 output is in kbar
 """

 f=open(outVASP,'r')
 word=[]
 sigma=np.zeros((3,3))
 i=0
 for line in f:
   if re.match('(.*)in kB(.*)', line):
     i+=1
 icount=i
 f.seek(0)
 i=0
 for line in f:
   if re.match('(.*)in kB(.*)', line):
     i+=1
     if (i==icount):
       word=line.split()
       sigma[0,0]=float(word[2])
       sigma[1,1]=float(word[3])
       sigma[2,2]=float(word[4])

       sigma[0,1]=float(word[5])
       sigma[1,2]=float(word[6])
       sigma[2,0]=float(word[7])

       sigma[1,0]=sigma[0,1]
       sigma[2,1]=sigma[1,2]
       sigma[0,2]=sigma[2,0]
       break
 f.close()
#patched version ... no conversion in kbar, already provided in kbar by VASP:
# sigma=sigma*1.e-3
 return sigma;

def get_sigma_from_vasp(outVASP, Patched=None):
 """
 Get sigma from patched output of VASP
 output is in kbar
 """
 sigma=np.zeros((3,3))
 word=[]
 with open(outVASP,'r') as f:
   m=mmap.mmap(f.fileno(),0,prot=mmap.PROT_READ)
   i=m.rfind(b"in kB")
   m.seek(i)
   line=m.readline()
   word=line.split()
   sigma[0,0]=float(word[2])
   sigma[1,1]=float(word[3])
   sigma[2,2]=float(word[4])

   sigma[0,1]=float(word[5])
   sigma[1,2]=float(word[6])
   sigma[2,0]=float(word[7])

   sigma[1,0]=sigma[0,1]
   sigma[2,1]=sigma[1,2]
   sigma[0,2]=sigma[2,0]
 f.close()
#patched version ... no conversion in kbar, already provided in kbar by VASP:
# sigma=sigma*1.e-3
 return sigma;


def get_volume_from_vasp_outcar_old(outVASP):
 f=open(outVASP,'r')
 mark=0
 itest=0
 abox = []
 bbox = []
 cbox = []
 for line in f:
   if re.match("(.*)direct lattice vectors(.*)", line):
     word=line.split()
     mark=1
     #print line
   if (mark>=1) and (mark<=3):
    if (itest==1):
     word=line.split()
     for i in range(3):
      abox.append(float(word[i]))
      #abox[i]=float(word[i])
    if (itest==2):
     word=line.split()
     for i in range(3):
      bbox.append(float(word[i]))
      #bbox[i]=float(word[i])
    if (itest==3):
     word=line.split()
     for i in range(3):
      cbox.append(float(word[i]))
      #cbox[i]=float(word[i])
    itest=itest+1
   marx=mark+1
   if (mark==4):
    break

 vect = crossProdNew(bbox,cbox)
 volume=dotProd(abox,vect)
 f.close()
 return volume;


def get_volume_from_vasp_outcar(outVASP, box=False):

 abox=[]
 bbox=[]
 cbox=[]

 with open(outVASP,'r') as f:
   m=mmap.mmap(f.fileno(),0,prot=mmap.PROT_READ)
   i=m.find(b"direct lattice vectors")
   m.seek(i+1)
   line=m.readline()

   line=m.readline()
   word=line.split()
   for i in range(3):
      abox.append(float(word[i]))

   line=m.readline()
   word=line.split()
   for i in range(3):
      bbox.append(float(word[i]))

   line=m.readline()
   word=line.split()
   for i in range(3):
      cbox.append(float(word[i]))

 if box==True:
   a= np.array(abox)
   b= np.array(bbox)
   c= np.array(cbox)
   box_vasp=np.vstack([ a, b, c])


 vect = crossProdNew(bbox,cbox)
 #this is gin volume ... 
 volume=dotProd(abox,vect)
 #this is bcc volume 
 #volume=2.0*dotProd(abox,vect)

 #debug_vol print("vol: ", abox, bbox, cbox, volume)
 f.close()
 if box==False:
    return volume;
 if box==True:
    return volume, box_vasp;


#------------------
# WARNINF THIS IS NOR CORRECT. replaced by new one
def crossProd(a,b):
# crros or vectorial product in n dimension
 dimension = len(a)
 c = []
 for i in range(dimension):
    c.append(0)
    for j in range(dimension):
      if j != i:
        for k in range(dimension):
          if k != i:
            if k > j:
              c[i] += a[j]*b[k]
            elif k < j:
              c[i] -= a[j]*b[k]
 return c

def crossProdNew(a,b):
# crros or vectorial product in n dimension
 if (  len(a) != len(b) ) :
   print('the vectors a and b do not have the same dimension. I m confused')
   exit(0)
 c = []
 c = np.cross(a,b)
 return c





def dotProd(a,b):
# dor or scalar product in n dimension
 dimension = len(a)
 c = 0.0
 for i in range(dimension):
  c  += a[i]*b[i]
 return c



def gen_gpfile_for_parabola_fit(a,b,c,prefix):
#Generate input script for GNUplot for fitting
  os.system("rm -f fit.log")
  fn = "%s_fit.gp"%prefix
  fout = open(fn,"w")
  fout.write('set term gif\n')
  fout.write('set output "%s_fit.gif"\n'%prefix)
  fout.write("f(x) = a*x**2 + b*x+ c\n")
  b = 0.0
  fout.write("a = %f; b = %f; c = %f\n"%(a,b,c))
  fout.write("FIT_LIMIT = 1e-9\n")
  if b == 0.0:
    fout.write("fit f(x) '%s_summary' via a,c\n"%prefix)
  else:
    fout.write("fit f(x) '%s_summary' via a,b,c\n"%prefix)
  fout.write('fx_title = sprintf("(%f)*x**2+(%f)*x+(%f)",a,b,c)\n')
  fout.write("plot '%s_summary' pt 5 t 'DATA',f(x) lc rgb 'blue' t fx_title"%prefix)
  fout.close()
  return fn

def get_fitted_param(a,b,c,fn):
#Extract fitted parameters from GNUplot output
  fin = file(fn,"r")
  while True:
    line = fin.readline()
    if line == "": sys.exit(123)
    if re.search("^Final set of parameters",line): break
  for i in range(3):
    line = fin.readline().split()
    #debug print line
    if   len(line)>0 and line[0] == "a": a = float(line[2])
    elif len(line)>0 and line[0] == "b": b = float(line[2])
    elif len(line)>0 and line[0] == "c": c = float(line[2])
  fin.close()
  #debug print 'in abc', a ,b ,c
  return a,b,c


def get_nionic_step_from_vasp_outcar(outVASP):
 f=open(outVASP,'r')
 nionic=0
 for line in f:
  if re.match('(.*) TOTAL-FORCE (.*)', line):
    nionic+=1
 f.close()
 return nionic;


#>--------------------------
def test_if_ndm_terminated_correctly(outVASP):
# return 0 is the OUTCAR is incomplete
# return 1 is the run is complete
 f=open(outVASP,'r')
 itest=1
 i=0
 for line in f:
  i=i+1
  last = line.split()
  if last[-1]=='NaN':
    itest=0
  if i > 1:
     break
 f.close()
 return itest;






def test_if_vasp_terminated_correctly(outVASP):
# return 0 is the OUTCAR is incomplete
# return 1 is the run is complete
 i=-100
 itest=0
 with open(outVASP,'r') as f:
   m=mmap.mmap(f.fileno(),0,prot=mmap.PROT_READ)
   i=m.rfind('Total CPU time used'.encode('utf-8'))
   if (i>0):
     itest=1
 f.close()
 return itest;
#<------------------------------




def print_last_forces_from_vasp(outVASP,file_to_write_forces, nat,nionic):


  fout=open(outVASP,'r')
  fforces=open(file_to_write_forces,'a')

  #print nat,nionic
  nstop=0
  nline=0
  narray=0
  ddrift=nionic*[3*[0.0]]
  ntest=0
  for line in fout:
    nline+=1

    if re.match('(.*) TOTAL-FORCE (.*)', line):
      nstop=1
      nstart=nline
      i=0
      ntest+=1
      #debug print nionic, ntest
    if (nionic==ntest):
      if (nstop==1):
        if ((nline>=(nstart+2)) and (nline<(nstart+2+nat))):
           #print line
           i+=1
           word=line.split()
           #print len(word)
           #print narray,i
           #if (i==1):
           fforces.write("%12.6f %12.6f %12.6f\n"%(float(word[3]), float(word[4]), float(word[5])))
           #debug(if you do not trust) print"%i %12.5f  %12.5f  %12.5f  %12.5f"%(nline, float(word[3]), float(word[4]), float(word[5]), ftemp)
        if (nline==(nstart+1+nat)):
           nstop==0
  fforces.close()
  return;


def print_structures_vasp(nat_per_box,nrepeat,boxfin,itype,rcell,cell):
  if np.amax(itype) > 2:
    print('WARNING: print_structure_vasp doesnt handle more than two types')
    exit(0)

  fout = open("POSCAR","w")
  fout.write("W\n")
  fout.write("1.0\n")
  fout.write(" %22.16f  %22.16f  %22.16f\n"%(boxfin[0,0], boxfin[0,1],boxfin[0,2]))
  fout.write(" %22.16f  %22.16f  %22.16f\n"%(boxfin[1,0], boxfin[1,1],boxfin[1,2]))
  fout.write(" %22.16f  %22.16f  %22.16f\n"%(boxfin[2,0], boxfin[2,1],boxfin[2,2]))

  if np.amax(itype) == 2:
   fout.write("%i %i\n"%(np.count_nonzero(itype == 1), np.count_nonzero(itype ==2)))
  if np.amax(itype) == 1:
   fout.write("%i\n"%(nat_per_box))
  fout.write("Direct\n")
  for i in range(nat_per_box):
    fout.write("%15.7f %15.7f %15.7f   \n"%(cell[i,0],cell[i,1],cell[i,2]))
  fout.close()
  return
