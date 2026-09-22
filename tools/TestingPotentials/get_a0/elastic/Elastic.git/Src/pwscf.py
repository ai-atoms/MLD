import globalv
import os
import re
import sys
import math
import numpy as np
import setup_param_pwscf

def setup_pwscf(based,indir,cname):

  os.chdir(based+'/'+indir)
  setup_param_pwscf.init()
  incar=setup_param_pwscf.incar
  kpoints=setup_param_pwscf.kpoints
  computer=globalv.cluster
  prefix=setup_param_pwscf.prefix

  if computer=='zelda':
     generate_jsub_pwscf_in_zelda(indir,cname)
  if computer=='occigen':
     generate_jsub_pwscf_in_occigen(indir,cname)
  if computer=='irene':
     generate_jsub_pwscf_in_irene(indir,cname)

  fj=open(prefix+'.scf.in','w')

  fi=open(incar,'r')
  for line in fi:
    fj.write(line)
  fi.close()

  fk=open(kpoints,'r')
  for line in fk:
    fj.write(line)
  fk.close()

  fp=open('POSCAR.pwscf','r')
  for line in fp:
    fj.write(line)
  fp.close()

  fj.close()
# This is for BSC
#  generate_jsub_bsc(indir,cname)
#This is for Occigen
#  generate_jsub_occigen(indir,cname)
#This is for Curie


  return;

def run_pwscf(outPWSCF):

 script=globalv.script
 setup_param_pwscf.init()
 computer=globalv.cluster


 os.system('pwd')
 fj=open('jsub','a')
 fj.write('mv fe.scf.out %s'%outPWSCF)
 fj.close()


 if computer=='zelda':
   run_pwscf_zelda(outPWSCF)

 if computer=='occigen':
   run_pwscf_occigen(outPWSCF)

 if computer=='irene':
   run_pwscf_irene(outPWSCF)



 if computer=='my_computer':
   run_ndm_my_computer(outPWSCF)
 return

def run_pwscf_zelda(outNDM):

#  if not os.path.exists(exerun):
#   print "Exe file %s doesn't exist. Put the correct path in setup_vasp"%str(exerun)

  os.system('pwd')
  os.system('qsub jsub')
  return

def run_pwscf_occigen(outNDM):

#  if not os.path.exists(exerun):
#   print "Exe file %s doesn't exist. Put the correct path in setup_vasp"%str(exerun)

  os.system('pwd')
  os.system('sbatch  jsub')
  return

def run_pwscf_irene(outNDM):

#  if not os.path.exists(exerun):
#   print "Exe file %s doesn't exist. Put the correct path in setup_vasp"%str(exerun)

  os.system('pwd')
  os.system('ccc_msub jsub')
  return




def generate_jsub_pwscf_in_occigen(indir,cname):

   setup_param_pwscf.init()
   script=globalv.script
   cal_elas=globalv.cal_elas
   PWX=setup_param_pwscf.PWX
   nodes=setup_param_pwscf.nodes
   nprocs=setup_param_pwscf.ntasks
   npools=setup_param_pwscf.npools
   prefix=setup_param_pwscf.prefix
   if (nprocs < 0):
     nprocs=int(nodes)*16

   f=open('jsub','w')
   f.write('#!/bin/sh\n')
   f.write('#SBATCH -J %s%s%s\n'%(prefix,indir,cname))
   f.write('#SBATCH --nodes=%s                # debug:  number of nodes \n'%nodes)
   f.write('#SBATCH --ntasks-per-node=16    # debug:  tasks /  nodes \n')
   f.write('#SBATCH --ntasks=%s              # the number of MPI Tasks \n'%nprocs)
   f.write('#SBATCH --threads-per-core=1     #  \n')
   f.write('#SBATCH --time  20:10:55         # Time limit: HH:MM:SS format\n')
   f.write('##SBATCH -p highmem              #  By default 2000 MB are reserved. More uncomment\n')
   f.write('#SBATCH --output  ofile.out   # stdout pf the job\n')
   f.write('#SBATCH --error   efile.err   # stderr of the job\n')
   f.write('\n')
   f.write('pwd\n')
   f.write('## SLURM_SUBMIT_DIR = the directory from which sbatch was invoked \n')
   f.write('## SLURM_JOB_ID  = job ID of the executing job \n')
   f.write('\n')
   f.write('module load espresso\n')
   f.write('module list\n')
   f.write('\n')
   f.write('suff=%s%s\n'%(indir,cname))
   f.write('export nprocs=${SLURM_NTASKS}\n')
   f.write('export npools=%s\n'%npools)
   f.write('\n')
   f.write('SUBMISSION=${SLURM_SUBMIT_DIR}\n')
   f.write('export OMP_NUM_THREADS=1\n')
   f.write('export PWX=%s\n'%PWX)
   f.write('SCRATCHDIR=$SCRATCHDIR\n')
   f.write('RUNDIR=${SCRATCHDIR}/%s/${suff}\n'%prefix)
   f.write('PWSCF_PSEUDO=$HOME/pwscf_work/pwscf_pseudo\n')
   f.write('rm -rf ${RUNDIR}\n')
   f.write('mkdir -p ${RUNDIR}\n')
   f.write('cd $RUNDIR\n')
   f.write('pwd\n')
  #f.write('cat > $LAUNCHDIR/fe.scf.in << ***\n')
   f.write("mpirun -np $nprocs $PWX  -nk $npools < $SUBMISSION/%s.scf.in > $SUBMISSION/%s.scf.out\n"%(prefix,prefix))
   f.write("cp $SUBMISSION/%s.scf.out $SUBMISSION/out.run\n"%(prefix))
   f.write("\n\n\n\n")
   if cal_elas=='qha_bulk':
      PHX=setup_param_pwscf.PHX
      ph_in=setup_param_pwscf.ph_in
      phrecover_in=setup_param_pwscf.phrecover_in

      Q2RX=setup_param_pwscf.Q2RX
      q2r_in=setup_param_pwscf.q2r_in

      MATDYNX=setup_param_pwscf.MATDYNX
      matdyn_in=setup_param_pwscf.matdyn_in

      nimages=setup_param_pwscf.nimages

      f.write('#PHonon part ...............\n')
      f.write('export nimages=%s\n'%nimages)
      f.write('export PHX=%s\n'%PHX)
      f.write('export Q2RX=%s\n'%Q2RX)
      f.write('export MATDYNX=%s\n'%MATDYNX)

      f.write("mpirun -np $nprocs $PHX  -ni $nimages -nk $npools < %s > $SUBMISSION/%s.ph.out\n"%(ph_in,prefix))
      newprocs=int(nprocs)/int(nimages)
      f.write("mpirun -np %i $PHX  -nk $npools < %s > $SUBMISSION/%s.ph_recover.out\n"%(newprocs,phrecover_in,prefix))
      f.write("%s  < %s > %s.q2r.out\n"%(Q2RX,q2r_in,prefix))
      f.write("mpirun -np 4 %s  < %s > $SUBMISSION/%s.matdyn.out\n"%(MATDYNX,matdyn_in,prefix))
      f.write("cp  matdyn.modes   $SUBMISSION/%s.matdyn.modes\n"%(prefix))

   f.close()
   return;


def generate_jsub_bsc(indir,cname):


   f=open('jsub','w')
   f.write('#!/bin/bash\n')
   f.write('#BSUB -J %s%s\n'%(indir,cname))
   f.write('##BSUB -cwd pathname # working directory. If not specified it the current workingdirectory\n')
   f.write('#BSUB -q prace       # debug:  queue max 64c/4n and 1h, MAX 1job/user prace: 4096\n')
   f.write('#BSUB -n 16          # The number of MPI tasks \n')
   f.write('#BSUB -W 00:55       # Time limit: HH:MM format\n')
   f.write('#BSUB -M 1800        # The maximum amount of memory necessary for each task (in MB).  By default 1800 MB are reserved\n')
   f.write('##BSUB -R "span[ptile=8]"  # If you want 4 process assigned / node. Any number between 1 and 16. \n')
   f.write('##                         # USE THIS OPTION WITH export OMP_NUM_THREADDS=4 (if you want threads)\n')
   f.write('#BSUB -oo file.o%J   # stdout pf the job\n')
   f.write('#BSUB -eo file.e%J   # stderr of the job\n')
   f.write('\n')
   f.write('\n')
   f.write('module purge\n')
   f.write('module load intel/13.0.1\n')
   f.write('module load impi/4.1.3.049\n')
   f.write('module load MKL/11.1.2\n')
   f.write('module load FFTW/3.3\n')
   f.write('\n')
   f.write('# Loadinh IntelMPI module load impi\n')
   f.write('#module load intel impi \n')
   f.write('\n')
   f.write('\n')
   f.write('module list\n')
   f.write('\n')
   f.write('# Loaging IBM MPI: module load poe\n')
   f.write('#module list\n')
   f.write('## module load VASP/5.3.5\n')
   f.write('suff=%s\n'%(indir))
   f.write('#export nprocs=${BRIDGE_MSUB_NPROC}\n')
   f.write('export nprocs=${LSB_DJOB_NUMPROC}\n')
   f.write('# there is also $LSB_MAX_NUM_PROCESSORS=64\n')
   f.write('#set > out.set\n')
   f.write('echo ${nprocs} \n')
   f.write('\n')
   f.write('SUBMISSION=${LS_SUBCWD}\n')
   f.write('#export OMP_NUM_THREADS=1\n')
   f.write('export EXE=${HOME}/Vasp/bin/vasp_kpoint_stand_ii.x\n')
   f.write('#export EXE=/apps/VASP/5.3.5/OPENMPI/bin/vasp.complex\n')
   f.write('export VASP_PSEUDO=${SUBMISSION}/POTCAR\n')
   f.write('export VASP_TMPDIR=${SCRATCHDIR}/TestPP/%s/${suff}\n'%(cname))
   f.write('\n')
   f.write('rm -rf   ${VASP_TMPDIR}\n')
   f.write('mkdir -p ${VASP_TMPDIR}\n')
   f.write('cd ${VASP_TMPDIR}\n')
   f.write('\n')
   f.write('pwd\n')
   f.write('\n')
   f.write('#SUBMISSION=${BRIDGE_MSUB_PWD}\n')
   f.write('#LSB_JOB_CWD: specifies the current directory for job execution\n')
   f.write('#LSB_HOSTS the list of hosts used to run the batch job \n')
   f.write('\n')
   f.write('\n')
   f.write('cp ${SUBMISSION}/INCAR   .\n')
   f.write('cp $VASP_PSEUDO               .\n')
   f.write('cp ${SUBMISSION}/POSCAR  .\n')
   f.write('cp ${SUBMISSION}/KPOINTS .\n')
   f.write('\n')
   f.write('outputfile=${SUBMISSION}/w${suff}.out\n')
   f.write('ulimit -c unlimited\n')
   f.write('#export MPI_GROUP_MAX=1024\n')
   f.write('\n')
   f.write('mpirun  -n $nprocs ${EXE}   |tee $outputfile \n')
   f.write('\n')
   f.write('cp OUTCAR  ${SUBMISSION}/\n')
   f.write('cp OSZICAR ${SUBMISSION}/\n')
   return;


def generate_jsub_pwscf_in_irene(indir,cname):

  setup_param_pwscf.init()
  script=globalv.script
  cal_elas=globalv.cal_elas
  PWX=setup_param_pwscf.PWX
  prefix=setup_param_pwscf.prefix
  nprocs=setup_param_pwscf.ntasks
  npools=setup_param_pwscf.npools
  prefix=setup_param_pwscf.prefix

  f=open('jsub','w')
  f.write('#!/bin/bash\n')
  f.write('##MSUB -q ivybridge \n')
  f.write('##MSUB -q broadwell \n')
  f.write('#MSUB -q standard \n')
  f.write('##MSUB -q large\n')
  f.write('#MSUB -r j%s%s%s\n'%(prefix,indir,cname))
  f.write('#MSUB -n %s    \n'%nprocs)
  f.write('#MSUB -T 86400   \n')
  f.write('#MSUB -o file.o%J \n')
  f.write('#MSUB -e file.e%J  \n')
  f.write('##MSUB -A den \n')
  f.write('#MSUB -A genden \n')
  f.write('##MSUB -A gen6973\n')
  f.write('\n')
  f.write('suff=%s%s\n'%(indir,cname))
  f.write('module load espresso/5.1.1\n')
  f.write('\n')
  f.write('export nprocs=${BRIDGE_MSUB_NPROC}\n')
  f.write('export npools=%s\n'%npools)
  f.write('export OMP_NUM_THREADS=1\n')
  f.write('export SUBMISSION=${BRIDGE_MSUB_PWD}\n')
  f.write('PWX=%s\n'%PWX)

  f.write('export PWSCF_PSEUDO=${CCCWORKDIR}/pwscf_pseudo\n')
  f.write('export RUNDIR=${CCCSCRATCHDIR}/%s/${suff}\n'%prefix)
  f.write('rm -rf ${RUNDIR}\n')
  f.write('mkdir -p ${RUNDIR}\n')
  f.write('cd ${RUNDIR}\n')
  f.write('pwd\n')
  f.write('#SCF part ...............\n')
  f.write('mpirun -np ${nprocs} $PWX  -nk ${npools}  < $SUBMISSION/%s.scf.in > $SUBMISSION/%s.scf.out\n'%(prefix,prefix))
  f.write("cp $SUBMISSION/%s.scf.out $SUBMISSION/out.run\n"%(prefix))
  f.write("\n\n\n\n")
  if cal_elas=='qha_bulk':
      PHX=setup_param_pwscf.PHX
      ph_in=setup_param_pwscf.ph_in
      phrecover_in=setup_param_pwscf.phrecover_in

      Q2RX=setup_param_pwscf.Q2RX
      q2r_in=setup_param_pwscf.q2r_in

      MATDYNX=setup_param_pwscf.MATDYNX
      matdyn_in=setup_param_pwscf.matdyn_in

      nimages=setup_param_pwscf.nimages

      f.write('#PHonon part ...............\n')
      f.write('export nimages=%s\n'%nimages)
      f.write('export PHX=%s\n'%PHX)
      f.write('export Q2RX=%s\n'%Q2RX)
      f.write('export MATDYNX=%s\n'%MATDYNX)

      f.write("mpirun -np $nprocs $PHX  -ni $nimages -nk $npools < %s > $SUBMISSION/%s.ph.out\n"%(ph_in,prefix))
      newprocs=int(nprocs)/int(nimages)
      f.write("mpirun -np %i $PHX  -nk $npools < %s > $SUBMISSION/%s.ph_recover.out\n"%(newprocs,phrecover_in,prefix))
      f.write("%s  < %s > %s.q2r.out\n"%(Q2RX,q2r_in,prefix))
      f.write("mpirun -np 4 %s  < %s > $SUBMISSION/%s.matdyn.out\n"%(MATDYNX,matdyn_in,prefix))
      f.write("cp  matdyn.modes   $SUBMISSION/%s.matdyn.modes\n"%(prefix))

  f.close()
  return;


def generate_jsub_pwscf_in_zelda(indir,cname):
#this jsub for zelda my friends ...

  setup_param_pwscf.init()
  script=globalv.script
  cal_elas=globalv.cal_elas
  PWX=setup_param_pwscf.PWX
  prefix=setup_param_pwscf.prefix
  nprocs=setup_param_pwscf.ntasks
  npools=setup_param_pwscf.npools
  prefix=setup_param_pwscf.prefix


  f=open('jsub','w')
  f.write('#!/bin/sh\n')
  f.write('#PBS -q prod_para\n')
  f.write('#PBS -l select=1:ncpus=%s\n'%nprocs)
  f.write('#PBS -l walltime=10:00:00\n')
  f.write('#PBS -N %s%s%s\n'%(prefix, indir,cname))
  f.write('#PBS -e error\n')
  f.write('#PBS -o output\n')
  f.write('\n')
  #f.write('nprocs=$PBS_NNODES\n')
  f.write('nprocs=%s\n'%nprocs)
  f.write('echo $PBS_NNODES\n')
  f.write('export pools=%s\n'%npools)
  f.write('export OMP_NUM_THREADS=1\n')
  f.write('cd $PBS_O_WORKDIR\n')
  f.write('pwd\n')
  f.write('suff=%s%s\n'%(indir,cname))

  f.write('PWX=/soft/QE/QE-5.4.0/espresso-5.4.0/bin/pw.x\n')
  f.write('PHX=/soft/QE/QE-5.4.0/espresso-5.4.0/bin/ph.x\n')
  f.write('Q2RX=/soft/QE/QE-5.4.0/espresso-5.4.0/bin/q2r.x\n')
  f.write('MATDYNX=/soft/QE/QE-5.4.0/espresso-5.4.0/bin/matdyn.x\n')

  f.write('SUBMISSION=$PBS_O_WORKDIR\n')
  f.write('SCRATCHDIR=/scratch/marinica/CalcPP\n')
  f.write('PWSCF_PSEUDO=$HOME/pwscf_work/pwscf_pseudo\n')
  f.write('RUNDIR=${SCRATCHDIR}/%s/${suff}\n'%prefix)
  f.write('rm -rf ${RUNDIR}\n')
  f.write('mkdir -p ${RUNDIR}\n')
  f.write('cd $RUNDIR\n')
  f.write('pwd\n')
  #f.write('cat > $LAUNCHDIR/fe.scf.in << ***\n')
  f.write("mpirun -np $nprocs $PWX  -nk $npools < $SUBMISSION/%s.scf.in > $SUBMISSION/%s.scf.out\n"%(prefix,prefix))
  f.write("cp $SUBMISSION/%s.scf.out $SUBMISSION/out.run\n"%(prefix))
  f.write("\n\n\n\n")
  if cal_elas=='qha_bulk':
      PHX=setup_param_pwscf.PHX
      ph_in=setup_param_pwscf.ph_in
      phrecover_in=setup_param_pwscf.phrecover_in

      Q2RX=setup_param_pwscf.Q2RX
      q2r_in=setup_param_pwscf.q2r_in

      MATDYNX=setup_param_pwscf.MATDYNX
      matdyn_in=setup_param_pwscf.matdyn_in

      nimages=setup_param_pwscf.nimages

      f.write('#PHonon part ...............\n')
      f.write('export nimages=%s\n'%nimages)
      f.write('export PHX=%s\n'%PHX)
      f.write('export Q2RX=%s\n'%Q2RX)
      f.write('export MATDYNX=%s\n'%MATDYNX)

      f.write("mpirun -np $nprocs $PHX  -ni $nimages -nk $npools < %s > $SUBMISSION/%s.ph.out\n"%(ph_in,prefix))
      newprocs=int(nprocs)/int(nimages)
      f.write("mpirun -np %i $PHX  -ni $nimages -nk $npools < %s > $SUBMISSION/%s.ph_recover.out\n"%(newprocs,phrecover_in,prefix))
      f.write("%s  < $SUBMISSION/%s > %s.q2r.out\n"%(q2r_in,prefix,prefix))
      f.write("mpirun -np 4 %s  < %s > $SUBMISSION/%s.matdyn.out\n"%(MATDYNX,matdyn_in,prefix))
      f.write("cp  matdyn.modes   $SUBMISSION/%s.matdyn.modes\n"%(prefix))


  f.close()
  return

def get_energy_from_pwscf_outcar(outPWSCF):

 if not os.path.exists(outPWSCF):
   print("%s  output of pwscf  doesn't exist. "%str(outPWSCF))
   os.system('pwd')
   exit(0)

 ry2ev=13.605698
 f=open(outPWSCF,'r')
 word=[]
 for line in f:
   if re.match("(.*)!    total energy (.*)", line):
     word=line.split()
     break
     #print line
 f.close()
 if len(word)==0:
   print('get_energy_from_pwscf_outcar: There is no energie output in the directory %s'%os.getcwd())
   exit(0)
 ene=float(word[4])*float(ry2ev)
 return ene;


def get_volume_from_pwscf_outcar(outPWSCF):
 b2a=0.5291772109
 f=open(outPWSCF,'r')
 for line in f:
   if re.match("(.*)unit-cell volume(.*)", line):
     word=line.split()
     vol=float(word[3])
     break
 f.close()
 if len(word)==0:
   print('get_volume_from_pwscf_outcar: There is no energie output in the directory %s'%os.getcwd())
   exit(0)

 volume=vol*b2a**3
 return volume;




def get_volume_old_from_pwscf_outcar(outPWSCF):
 b2a=0.5291772109
 f=open(outPWSCF,'r')
 mark=0
 itest=0
 abox = []
 bbox = []
 cbox = []
 for line in f:

   if re.match("(.*)lattice parameter(.*)", line):
     word=line.split()
     if len(word)==0:
       print('get_volume_from_pwscf_outcar: There is no alat in  output in the directory %s'%os.getcwd())
       exit(0)
     alat=float(word[4])*b2a
   if re.match("(.*) crystal axes(.*)", line):
     mark=1
   if (mark>1) and (mark<=4):

    if (mark==2):
     word=line.split()
     for i in range(3,6):
      abox.append(float(word[i]))
      #abox[i]=float(word[i])
     #print abox
    if (mark==3):
     word=line.split()
     for i in range(3,6):
      bbox.append(float(word[i]))
      #bbox[i]=float(word[i])
     #print bbox
    if (mark==4):
     word=line.split()
     for i in range(3,6):
      cbox.append(float(word[i]))
      #cbox[i]=float(word[i])
     #print cbox
    itest=itest+1
   if (mark>0):
     mark=mark+1
     if (mark==5):
      break
 abox=np.array(abox)*alat
 bbox=np.array(bbox)*alat
 cbox=np.array(cbox)*alat
 vect = crossProdNew(bbox,cbox)
 volume=dotProd(abox,vect)
 f.close()
 return volume;


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

def read_eigenvalues_from_pwscf(file_eigenvalues):
  icount=0
  icountq=0
  omega=[]
  word=[]
  feigen=open(file_eigenvalues,'r')
  for line in  feigen:
     if re.match("(.*) freq (.*)",line):
         word=line.split()
         icount=icount+1
         #print word[4]
         omega.append(float(word[4]))
     if re.match("(.*) q = (.*)",line):
         icountq=icountq+1

  omega=np.array(omega)
  #print icount , omega.size
  if len(word)==0:
   print('read_eigenvalues_from_pwscf: There is no energie output in the directory %s'%os.getcwd())
   exit(0)
  return icountq, omega


def get_natom_from_pwscf_outcar(outPWSCF):
 if not os.path.exists(outPWSCF):
   print("%s  output of pwscf  doesn't exist. "%str(outPWSCF))
   os.system('pwd')
   exit(0)

 f=open(outPWSCF,'r')
 for line in f:
   if re.match("(.*)number of atoms(.*)", line):
     word=line.split()
     nat=float(word[4])
     break
 f.close()
 if len(word)==0:
   print('get_natom_from_pwscf_outcar: There is no energie output in the directory %s'%os.getcwd())
   exit(0)
 return nat;


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
  fin = open(fn,"r")
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

def print_structures_pwscf(nat_per_box,nrepeat,boxfin,itype,rcell,cell):
  if np.amax(itype) > 1:
    print('WARNING: print_structure_pwscf doent handle more than one type')
    exit(0)
  fout = open("POSCAR.pwscf","w")
  fout.write("CELL_PARAMETERS { angstrom } \n")
  fout.write(" %22.16f  %22.16f  %22.16f\n"%(boxfin[0,0], boxfin[0,1],boxfin[0,2]))
  fout.write(" %22.16f  %22.16f  %22.16f\n"%(boxfin[1,0], boxfin[1,1],boxfin[1,2]))
  fout.write(" %22.16f  %22.16f  %22.16f\n"%(boxfin[2,0], boxfin[2,1],boxfin[2,2]))
  #fout.write("%i\n"%(nat_per_box))
  fout.write("ATOMIC_POSITIONS { crystal }\n")
  for i in range(nat_per_box):
    fout.write("Fe   %15.7f %15.7f %15.7f   \n"%(cell[i,0],cell[i,1],cell[i,2]))
  fout.close()
  return
