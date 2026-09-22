import os
import re
import sys
import math
import globalv
import setup_param_lammps
import numpy as np
from numpy.linalg import norm
from tools import is_upper_triangular,  right_hand_basis
# setup_lammps    ... only for minimizations ...
# setup_ndm_ph    ... active for phondy calculation. Called from qha_volume.py

def setup_lammps(based,indir,cname):

  setup_param_lammps.init()
  os.chdir(based+'/'+indir)
  eampot=setup_param_lammps.eampot
  dinfile=setup_param_lammps.dinfile
  computer=globalv.cluster
  dynamics=globalv.dynamics
  gap=1
  os.system("ln -s %s  pot.fs"%(eampot))
  os.system("cp -p %s  cube.in"%(dinfile))
  if dynamics=='neb':
    dinfile_neb=setup_param_lammps.dinfile_neb
    os.system("cp -p %s cube_neb.in"%(dinfile_neb))

  if gap == 1:
      pass
     #os.system("ln -s   /home/marinica/TESTS/GAP_Potentials/6_W/gp.xml  gp.xml")
     #os.system("ln -s   /home/marinica/TESTS/GAP_Potentials/6_W/gp.xml.sparseX.GAP_2013_6_24_60_12_58_8_3271  gp.xml.sparseX.GAP_2013_6_24_60_12_58_8_3271")

  if computer=='zelda':
    generate_jsub_lammps_zelda(indir,cname)
  if computer=='occigen':
    generate_jsub_lammps_occigen(indir,cname)
  if computer=='irene':
    generate_jsub_lammps_occigen(indir,cname)



def run_lammps(based, dirrun, outNDM):
 script=globalv.script
 setup_param_lammps.init()
 computer=globalv.cluster
 exerun=setup_param_lammps.EXELAMMPS



 if computer=='zelda':
   run_lammps_zelda(outNDM)

 if computer=='irene':
   run_lammps_irene(outNDM)

 if computer=='occigen':
   run_lammps_occigen(outNDM)

 if computer=='my_computer':
   run_lammps_my_computer(outNDM,dirrun)
 return



def run_lammps_my_computer(outNDM,dirrun):
  script=globalv.script
  cal_elas=globalv.cal_elas
  setup_param_lammps.init()
  exerun=setup_param_lammps.EXELAMMPS
  ntasks=globalv.ntasks
  dynamics=globalv.dynamics
  #debug print 'r_my_computer', exerun
  if dynamics=='min':
     os.system("mpirun -np %s %s < cube.in > %s "%(ntasks, exerun,  outNDM))
  elif dynamics=='neb':
     if dirrun[1:5]=='bulk':
       os.system("mpirun -np %s %s < cube.in > %s "%(ntasks, exerun,  outNDM))
     else:
       input_neb_incar=globalv.input_neb_incar
       os.system("cp cube_deb.lmp   cube.lmp")
       os.system("cp cube_deb.lmp   cube_deb_ini.lmp")
       os.system("mpirun -np 1 %s < cube.in > /dev/null"%(exerun))
       os.system("cp relaxed.data cube_deb.lmp")
       os.system("cp cube_fin.lmp   cube.lmp")
       os.system("cp cube_fin.lmp   cube_fin_ini.lmp")
       os.system("mpirun -np 1 %s < cube.in > /dev/null"%(exerun))
       os.system("cp relaxed.data cube_fin.lmp")
       nat_per_box, boxini, itype, cell, rcell = read_lammps_data('./',name_of_input_file='relaxed.data')
       fn=open('cube_fin.disp','w')
       fn.write('%i\n'%nat_per_box)
       for i in range(nat_per_box):
            fn.write(' %i %.15f %.15f %.15f\n'%(i+1, rcell[i,0],rcell[i,1],rcell[i,2]))
       fn.close()
       os.system("mpirun -np %s %s  -partition %sx1 -i  cube_neb.in > %s"%(ntasks, exerun, ntasks,   outNDM))
  else:
     print('Wrong dynamics 70')
     exit(0)

  if cal_elas=='qha_bulk':
    dinfile_cg=setup_param_ndm.dinfile_cg
    dinfile_ph=setup_param_ndm.dinfile_ph
    exerun_ph=setup_param_ndm.EXENDM_PH
    sauve2gin=setup_param_ndm.SAUVE2GIN
    relaxation=globalv.relaxation
    os.system(" cp -p %s cube.din"%(dinfile_cg))
    os.system("%s > out.run_cg"%exerun)
    if relaxation=='no':
      os.system('mv out.run_cg  out.run')
    elif relaxation=='yes':
      dinfile_cg2=setup_param_ndm.dinfile_cg2
      dinfile_tr=setup_param_ndm.dinfile_tr
      imm = get_imm_from_din('cube.din')
      os.system("%s cube.cout cube.gin <<< %s >/dev/null"%(sauve2gin, imm))

      os.system(" cp -p %s cube.din"%(dinfile_tr))
      os.system("%s > out.run_tr"%exerun)
      imm = get_imm_from_din('cube.din')
      os.system("%s cube.cout cube.gin <<< %s >/dev/null"%(sauve2gin, imm))


      os.system(" cp -p %s cube.din"%(dinfile_cg2))
      os.system("%s > out.run"%exerun)
      imm = get_imm_from_din('cube.din')
      os.system("%s cube.cout cube.gin <<< %s >/dev/null"%(sauve2gin, imm))
    else:
      print('wrong value for relaxation variable. Only yes and no allowed')
      exit(0)

    os.system(" cp -p %s cube.din"%(dinfile_ph))
    os.system(" mpirun -np 12  %s > %s "%(exerun_ph, 'out.run_ph'))

  return


def run_lammps_zelda(outNDM):

#  if not os.path.exists(exerun):
#   print "Exe file %s doesn't exist. Put the correct path in setup_vasp"%str(exerun)

  os.system('pwd')
  os.system('qsub jsub')
  return

def run_lammps_occigen(outNDM):

#  if not os.path.exists(exerun):
#   print "Exe file %s doesn't exist. Put the correct path in setup_vasp"%str(exerun)

  os.system('pwd')
  os.system('sbatch  jsub')
  return


def generate_jsub_lammps_occigen(indir,cname):
   setup_param_lammps.init()
   script=globalv.script
   EXELAMMPS=setup_param_lammps.EXELAMMPS
   potcar=setup_param_lammps.eampot
   incar=setup_param_lammps.dinfile
   nodes=setup_param_lammps.nodes
   nprocs=setup_param_lammps.ntasks
   nkpar=setup_param_lammps.nkpar
   type_node=setup_param_lammps.type_node
   prefix=globalv.prefix
   if (nprocs < 0):
     if (type_node=='HSW24'):
             task_per_node=24
             nprocs=int(nodes)*task_per_node
     elif (type_node=='BDW28'):
             task_per_node=28
             nprocs=int(nodes)*task_per_node
     else:
        print("Unknown type of node. For occigen it should be used only:")
        print("        HSW24 or BRW28")
        exit(0)

   f=open('jsub','w')
   f.write('#!/bin/sh\n')
   f.write('#SBATCH -J %s%s\n'%(indir,cname))
   f.write('#SBATCH --nodes=%s                # debug:  number of nodes \n'%nodes)
   f.write('#SBATCH --constraint=%s          # debug:  type of nodes \n'%type_node)
   f.write('#SBATCH --ntasks-per-node=%s      # debug:  tasks /  nodes \n'%task_per_node)
   f.write('#SBATCH --ntasks=%s               # the number of MPI Tasks \n'%nprocs)
   f.write('#SBATCH --threads-per-core=1      #  \n')
   f.write('#SBATCH --time  00:10:55         # Time limit: HH:MM:SS format\n')
   f.write('##SBATCH -p highmem              #  By default 2000 MB are reserved. More uncomment\n')
   f.write('#SBATCH --output  ofile.out   # stdout pf the job\n')
   f.write('#SBATCH --error   efile.err   # stderr of the job\n')
   f.write('\n')
   f.write('\n')
   f.write('module purge\n')
   f.write('module load intel\n')
   f.write('module load intelmpi/2017.0.098\n')
   f.write('module load python\n')
   f.write('##---------------------\n')
   f.write('pwd\n')
   f.write('## SLURM_SUBMIT_DIR = the directory from which sbatch was invoked \n')
   f.write('## SLURM_JOB_ID  = job ID of the executing job \n')
   f.write('\n')
   f.write('module list\n')
   f.write('\n')
   f.write('suff=%s%s\n'%(indir,cname))
   f.write('export nprocs=${SLURM_NTASKS}\n')
   f.write('echo ${nprocs} \n')
   f.write('\n')
   f.write('SUBMISSION=${SLURM_SUBMIT_DIR}\n')
   f.write('#export OMP_NUM_THREADS=1\n')
   f.write('export EXE=%s\n'%EXELAMMPS)
   f.write('export VASP_TMPDIR=${SCRATCHDIR}/TestPP/%s/${suff}\n'%(cname))
   f.write('\n')
   f.write('rm -rf   ${VASP_TMPDIR}\n')
   f.write('mkdir -p ${VASP_TMPDIR}\n')
   f.write('cd ${VASP_TMPDIR}\n')
   f.write('\n')
   f.write('pwd\n')
   f.write('\n')
   f.write('\n')
   f.write('\n')
   f.write('cp %s          cube.in\n'%incar)
   f.write('ln -s %s       pot.fs \n'%potcar)
   f.write('ln -s %s       pos.lmp \n'%potcar)
   f.write('cp ${SUBMISSION}/POSCAR  .\n')
   f.write('\n')
   f.write('outputfile=${SUBMISSION}/w${suff}.out\n')
   f.write('ulimit -s unlimited\n')
   f.write('#export MPI_GROUP_MAX=1024\n')
   f.write('\n')
   f.write('srun  --mpi=pmi2 -K1 --resv-ports -n ${SLURM_NTASKS}  ${EXE} < cube.in  |tee $outputfile \n')
   f.write('\n')
   f.write('cp ${outputfile}  ${SUBMISSION}/out.run\n')
   return;

def generate_jsub_lammps_zelda(indir,cname):
#this jsub for zelda my friends ...
  script=globalv.script
  cal_elas=globalv.cal_elas
  setup_param_lammps.init()
  exerun=setup_param_lammps.EXELAMMPS

  potcar=setup_param_lammps.eampot
  incar=setup_param_lammps.dinfile
  nodes=setup_param_lammps.nodes
  nprocs=setup_param_lammps.ntasks
  nkpar=setup_param_lammps.nkpar
  dynamics=setup_param_lammps.dynamics
  if (nodes > 1):
    print("lammps.py: on Zelda nodes cannot be higher than 1. Put 1!")
    exit(0)
  if (nprocs > 24):
    print("lammps.py: on Zelda nprocs cannot be higher than 24.")
    exit(0)


  f=open('jsub','w')
  f.write('#!/bin/sh\n')


  if nprocs==1:
    f.write('#PBS -q prod\n')
  if nprocs  > 1:
    f.write('#PBS -q prod_para\n')
  if dynamics=='neb':
    f.write('#PBS -l select=%s:ncpus=%s\n'%(nodes,'12'))
  else:
    f.write('#PBS -l select=%s:ncpus=%s\n'%(nodes,nprocs))
  f.write('#PBS -l walltime=10:00:00\n')
  f.write('#PBS -N %s%s\n'%(indir,cname))
  f.write('#PBS -e error\n')
  f.write('#PBS -o output\n')
  f.write('\n')
  f.write('sleep 1\n')
  f.write('cd $PBS_O_WORKDIR\n')
  f.write('pwd\n')
  f.write('\n')
  f.write('cd $PBS_O_WORKDIR\n')
  f.write('\n')
  f.write('module load lammps\n')
  f.write('module load python2\n')
  f.write('export LD_LIBRARY_PATH=/soft/Conda/miniconda2/lib/:${LD_LIBRARY_PATH}\n')
  f.write('module list')

  f.write('\n')
  f.write('#phonons calculation ...\n')
  mpirunexe='/soft/openmpi-2.0.0/bin/mpirun'
  if script=='elastic.py':
    if (dynamics=='neb' and indir[1:5] != 'bulk' ):
       f.write("%s -np %s %s  -partition %sx1 -i  cube_neb.in > out.run \n "%(mpirunexe,  nprocs,  exerun, nprocs))
    else:
       f.write('%s -np %s %s < cube.in > out.run\n'%(mpirunexe,nprocs, exerun))

    f.write('sleep 5\n')
  if cal_elas=='qha_bulk':
    print("lammps.py: not yet implementation for qha and lammps mode.")
    exit(0)
    exerun_ph=setup_param_ndm.EXENDM_PH
    exerun=setup_param_ndm.EXENDM
    dinfile_cg=setup_param_ndm.dinfile_cg
    relaxation=globalv.relaxation

    f.write(" cp -p %s cube.din\n"%(dinfile_cg))
    f.write("%s > out.run_cg\n"%exerun)
    if relaxation=='no':
     f.write("cp -pf out.run_cg out.run\n")
    elif relaxation=='yes':
     dinfile_cg2=setup_param_ndm.dinfile_cg2
     dinfile_tr=setup_param_ndm.dinfile_tr
     ELASTIC=setup_param_ndm.ELASTIC
     sauve2gin=setup_param_ndm.SAUVE2GIN

     f.write("python %s/Src/get_imm.py cube.din > imm.tmp \n"%ELASTIC)
     f.write("%s cube.cout cube.gin < imm.tmp >/dev/null\n"%sauve2gin)

     f.write(" cp -p %s cube.din\n"%(dinfile_tr))
     f.write("%s > out.run_tr\n"%exerun)
     f.write("python %s/Src/get_imm.py cube.din > imm.tmp\n"%ELASTIC)
     f.write("%s cube.cout cube.gin < imm.tmp >/dev/null\n"%sauve2gin)


     f.write(" cp -p %s cube.din\n"%(dinfile_cg2))
     f.write("%s > out.run\n"%exerun)
     f.write("python %s/Src/get_imm.py cube.din > imm.tmp\n"%ELASTIC)
     f.write("%s cube.cout cube.gin < imm.tmp >/dev/null\n"%sauve2gin)

     f.write(" cp -p %s cube.din\n"%(dinfile_tr))
     f.write("%s > out.run_tr\n"%exerun)
     f.write("python %s/Src/get_imm.py cube.din > imm.tmp\n"%ELASTIC)
     f.write("%s cube.cout cube.gin < imm.tmp >/dev/null\n"%sauve2gin)


     f.write(" cp -p %s cube.din\n"%(dinfile_cg2))
     f.write("%s > out.run\n"%exerun)
     f.write("python %s/Src/get_imm.py cube.din > imm.tmp\n"%ELASTIC)
     f.write("%s cube.cout cube.gin < imm.tmp >/dev/null\n"%sauve2gin)

     f.write(" cp -p %s cube.din\n"%(dinfile_tr))
     f.write("%s > out.run_tr\n"%exerun)
     f.write("python %s/Src/get_imm.py cube.din > imm.tmp\n"%ELASTIC)
     f.write("%s cube.cout cube.gin < imm.tmp >/dev/null\n"%sauve2gin)


     f.write(" cp -p %s cube.din\n"%(dinfile_cg2))
     f.write("%s > out.run\n"%exerun)
     f.write("python %s/Src/get_imm.py cube.din > imm.tmp\n"%ELASTIC)
     f.write("%s cube.cout cube.gin < imm.tmp >/dev/null\n"%sauve2gin)

    else:
      print('wrong value for relaxation variable. Only yes and no allowed')
      exit(0)


  f.close()
  return

#------neb tools are here.
def return_last_line (outNDM):


  with open(outNDM,'r') as fneb:
    lines=fneb.read().splitlines()
    last_line=lines[-1]

  return last_line;

def get_climb_replica_from_lammps(outNDM, Patched=None):
 f=open(outNDM,'r')
 word=[]

 if Patched is None:
   print('get_climb_replica_from_lammps: The output %s should comes from NEB'%os.getcwd())
   exit(0)

#patched version ...

 if Patched is not None:
   f.seek(0)
   for line in f:
     if re.match('(.*)Climbing replica(.*)', line):
       word=line.split()
       no_replica=float(word[-1])
       break

   f.close()

   try:
       no_replica
   except NameError:
     print("get_climb_replica_from_lammps: Probably, the output of LAMMPS not comes from NEB")
     exit(0)

 return float(no_replica);



def get_all_from_neb_lammps(outNDM):
  rcoord=[]
  energy=[]
  line=return_last_line(outNDM)
  word=line.split()
  itest=0
  ilast=0
  for i in range(len(word)):
     if float(word[i])==0:
       itest=1
     if (itest >= 1):
       if (itest-1)%2 == 1:
           energy.append(word[i])
           if ilast==1:
             break
       if (itest-1)%2 == 0:
           rcoord.append(word[i])
           if float(word[i]) == 1:
             ilast=1
       itest=itest+1

  xp_vol=np.array(rcoord).astype(np.float)
  y_atomic=np.array(energy).astype(np.float)

  #debug plot_data (xp_vol, y_atomic)

  ymin=np.amin(y_atomic)
  min_loc=np.argmin(y_atomic)

  ymax=np.amax(y_atomic)
  max_loc=np.argmax(y_atomic)

  if (min_loc != 0):
    print ("NEB WARNING the 0 image is not the lowest energy along the path")

  isaddle=get_climb_replica_from_lammps('out.run', Patched='Yes')
  file_saddle='screen.%i'%isaddle
  nat = get_natom_from_lammps(file_saddle,Patched=True)
  energy = get_energy_from_lammps(file_saddle,Patched=True)
  volume = get_volume_from_lammps(file_saddle,Patched=True)
  sigma = get_sigma_from_lammps(file_saddle,Patched=True)

  return nat, energy, volume, sigma;
#net tools up here.

def get_natom_from_lammps(outNDM, Patched=None):
 f=open(outNDM,'r')
 word=[]
 i=0
 for line in f:
   if re.match('(.*) atoms(.*)', line):
     word=line.split()
     if word[-1] == 'atoms':
         natoms=int(word[-2])
         break

 if len(word)==0:
   print('get_natom_from_lammps: There is no nat output in the directory %s'%os.getcwd())
   exit(0)

#patched version ...

 if Patched is not None:
   f.seek(0)
   for line in f:
     if re.match('(.*)Number_of_atoms(.*)', line):
       word=line.split()
       p_natoms=float(word[-1])
       break

   f.close()

   try:
       p_natoms
   except NameError:
     print("get_natom_from_lammps: error in number of atoms. Probably, the output of LAMMPS is not patched")
     exit(0)

   if natoms != p_natoms :
    print('get_natom_from_lammps: incorrect number of atoms %s %s %s'%(natoms, p_atoms, os.getcwd()))
    exit(0)



 return float(natoms);



def get_energy_from_lammps(outNDM, Patched=None):
 #DEBUG print('------------------')
 #DEBUG flw = os.getcwd() 
 #DEBUG print(flw)

 f=open(outNDM,'r')
 word=[]
 i=0
 
 if Patched is None:
   for line in f:
     if i==1:
       word=line.split()
       p_energy=float(word[-1])
       break
     if re.match('(.*) Energy initial(.*)', line):
       word=line.split()
       i=i+1

   if len(word)==0:
     print('get_energy_from_ndm: There is no energy output in the directory %s'%os.getcwd())
     exit(0)


#patched version ...

 if Patched is not None:
   f.seek(0)
   for line in f:
     if re.match('(.*)Energy_box(.*)', line):
       word=line.split()
       p_energy=float(word[-1])
       break

   f.close()

   try:
       p_energy
   except NameError:
     print("get_energy_from_lammps: error in energy. Probably, the output of LAMMPS is not patched")
     exit(0)
#   if math.fabs(energy - p_energy) > 1.e-5 :
#    print 'get_natom_from_lammps: incorrect or inconsistency in energy %s %s %s'%(energy, p_energy, os.getcwd())
#    exit(0)

 return p_energy;


def get_volume_from_lammps(outNDM, Patched=None):
 f=open(outNDM,'r')
 word=[]

 if Patched is None:
   print('get_volume_from_lammps: The output %s should pe patched with Volume keyword'%os.getcwd())
   exit(0)

#patched version ...

 if Patched is not None:
   f.seek(0)
   for line in f:
     if re.match('(.*)Volume(.*)', line):
       word=line.split()
       volume=float(word[-1])
       break

   f.close()

   try:
       volume
   except NameError:
     print("get_volume_from_lammps: error in volume. Probably, the output of LAMMPS is not patched")
     exit(0)

 return float(volume);


def get_sigma_from_lammps(outNDM, Patched=None):
 """
 Get sigma from patched output of LAMMPS
 patched means that the component of the stress are labeled by sig_xx, sig_xy, sig_xz, sig_yz, sig_zz
 output is in kbar
 """

 f=open(outNDM,'r')
 word=[]
 i=0
 sigma=np.zeros((3,3))
 i=0
 for line in f:
   if re.match('(.*)sig_xx(.*)', line):
     word=line.split()
     sigma[0,0]=float(word[-1])
     i+=1
   if re.match('(.*)sig_yy(.*)', line):
     word=line.split()
     sigma[1,1]=float(word[-1])
     i+=1
   if re.match('(.*)sig_xy(.*)', line):
     word=line.split()
     sigma[0,1]=float(word[-1])
     sigma[1,0]=float(word[-1])
     i+=1
   if re.match('(.*)sig_xz(.*)', line):
     word=line.split()
     sigma[0,2]=float(word[-1])
     sigma[2,0]=float(word[-1])
     i+=1
   if re.match('(.*)sig_yz(.*)', line):
     word=line.split()
     sigma[1,2]=float(word[-1])
     sigma[2,1]=float(word[-1])
     i+=1
   if re.match('(.*)sig_zz(.*)', line):
     word=line.split()
     sigma[2,2]=float(word[-1])
     i+=1


 if i!=6:
   print("get_sigma_from_lammps: error in sigma. Probably, the output of LAMMPS is not patched")
   exit(0)

#patched version ...convert in kbar:
 sigma=sigma*1.e-3
 return sigma;


def convert_cell_from_gin_to_lammps(box_cell):
    """
    Convert a parallelpiped (forming right hand basis)
    to lower triangular matrix LAMMPS can accept. This
    function transposes cell matrix so the bases are column vectors
    """
    cell = np.matrix.transpose(box_cell)

    if not is_upper_triangular(cell):
        # rotate bases into triangular matrix
        tri_mat = np.zeros((3, 3))
        A = np.array([cell[0, 0], cell[1,0], cell[2,0]])
        B = np.array([cell[0, 1], cell[1,1], cell[2,1]])
        C = np.array([cell[0, 2], cell[1,2], cell[2,2]])
        #old version
        #A = np.array(cell[:, 0])
        #B = np.array(cell[:, 1])
        #C = np.array(cell[:, 2])

        #print 'A', A, np.shape(A)
        #print 'B', B, np.shape(B)
        #print 'C', C, np.shape(C)
        if not right_hand_basis(A,B,C):
             print("WARNING: your reper is not right handed.")
             print("WARNING: This is a critical issue. The LAMMPS results are wrong !!!!!")


        tri_mat[0, 0] = norm(A)
        Ahat = A / norm(A)
        AxBhat = np.cross(A, B) / norm(np.cross(A, B))
        tri_mat[0, 1] = np.dot(B, Ahat)
        tri_mat[1, 1] = norm(np.cross(Ahat, B))
        tri_mat[0, 2] = np.dot(C, Ahat)
        tri_mat[1, 2] = np.dot(C, np.cross(AxBhat, Ahat))
        tri_mat[2, 2] = norm(np.dot(C, AxBhat))

        # create and save the transformation for coordinates
        volume = np.linalg.det(box_cell)
        trans = np.array([np.cross(B, C), np.cross(C, A), np.cross(A, B)])
        trans = trans / volume
        coord_transform = np.dot(tri_mat , trans)

        return tri_mat, coord_transform
    else:
        return cell, None

def test_if_lammps_terminated_correctly(outLAMMPS):
# return 0 is the OUTCAR is incomplete
# return 1 is the run is complete
 f=open(outLAMMPS,'r')
 itest=0
 for line in f:
  if re.match('(.*)Energy_box(.*)', line):
    itest=1
 f.close()
 return itest;



def read_lammps_data(dirSTR,name_of_input_file=None):

 if name_of_input_file is None:
   sgin= dirSTR+'/'+ 'structure.lmp'
 else:
   sgin= dirSTR+'/'+ name_of_input_file

 if not os.path.exists(sgin):
    print("lammps file %s doesn't exist. Put the correct path"%str(sgin))
    exit(0)

 fgin=open(sgin,'r')
 lines_fgin=fgin.read().splitlines()
 icount=0
 itest=99999999999999999999999

 for i in range(len(lines_fgin)):
   if not lines_fgin[i]=='':
     word=lines_fgin[i].split()
     if i==2:
       if word[-1]=='atoms':
          nat_per_box=int(word[0])
          rcell=np.matrix(np.zeros((nat_per_box,3)))
          cell=np.matrix(np.zeros((nat_per_box,3)))
          itype=np.array(np.zeros(nat_per_box))
       else:
          print('error in line 3 of LAMMPS')
          exit(0)

     if i==3:
       if word[-1]=='types':
         ntype=int(word[0])
       else:
         print('error in line 4 of LAMMPS')
         exit(0)


     if i==5:
        xlo=float(word[0])
        xhi=float(word[1])
        if not(word[-1]=='xhi'):
          print('not legal format for xhi')
          exit(0)

     if i==6:
        ylo=float(word[0])
        yhi=float(word[1])
        if not(word[-1]=='yhi'):
          print('not legal format for yhi')
          exit(0)

     if i==7:
        zlo=float(word[0])
        zhi=float(word[1])
        if not(word[-1]=='zhi'):
          print('not legal format for zhi')
          exit(0)
        if not(word[-2]=='zlo'):
          print('not legal format for zlo')
          exit(0)





     if (i==8):
        xy=float(word[0])
        xz=float(word[1])
        yz=float(word[2])
        if not(word[-1]=='yz'):
          print('not legal format for yz')
          exit(0)
        if (xlo > 0.0) or (ylo > 0.0) or (zlo > 0.0):
           print('format unknown for xlo , ylo or zlo > 0')
           exit(0)
        #lammps_cell[0, 0] = xhi
        #lammps_cell[1, 1] = yhi
        #lammps_cell[2, 2] = zhi
        #lammps_cell[0, 1] = xy
        #lammps_cell[0, 2] = xz
        #lammps_cell[1, 2] = yz
        #debug print xlo,ylo,zlo,xhi,yhi,zhi,xy,xz,yz
        a,b,c,alpha,beta,gamma=convert_xyzlo_into_abc(xlo,ylo,zlo,xhi,yhi,zhi,xy,xz,yz)
        boxini= convert_abc_into_box_cell(a,b,c,alpha,beta,gamma)
        #debug print boxini
     if i==10:
       if word[0]=='Masses':
          itest=16
       elif word[0]=='Atoms':
          itest=11
          #debug print itest, lines_fgin[12]
       else:
          print('error in LAMMPS file at line 10')
          exit(0)
     if (i>=itest) and (i <= (itest+nat_per_box)):
         if word==[]:
              itest=itest+1
              break
         iatom=int(word[0])
         itype[iatom-1]=int(word[1])
         rcell[iatom-1,0]=float(word[2])
         rcell[iatom-1,1]=float(word[3])
         rcell[iatom-1,2]=float(word[4])



 fgin.close()
 boxini_inv=np.linalg.inv(boxini)
 cell = np.matmul(boxini_inv,rcell.T).T

 return nat_per_box, boxini, itype, cell, rcell;



def  convert_xyzlo_into_abc(xlo,ylo,zlo,xhi,yhi,zhi,xy,xz,yz):
# This function converts xlo,ylo... lammps into
# a, b, c, alpha, beta, gamma format.


 x_t = xhi-xlo
 y_t = yhi-ylo
 z_t = zhi-zlo

 a = x_t
 b = np.sqrt(y_t**2+xy**2)
 c = np.sqrt(z_t**2 + xz**2 + yz**2)
 alpha = np.arccos( (xy*xz+y_t*yz)/(b*c))
 beta  = np.arccos(xz/c)
 gamma=np.arccos(xy/b)

 return a,b,c,alpha,beta,gamma;



def  convert_abc_into_box_cell(a_size,b_size,c_size,alpha,beta,gamma):

#  This function converts conventional LAMMPS vectors
#   defined by a b c alpha beta gamma,
#   into a lower triangular matrix:
#         |  H(1,1)    0       0     |
#    H =  |  H(2,1)  H(2,2)    0     |
#         |  H(3,1)  H(3,2)  H(3,3)  |

 box_cell=np.zeros((3,3))
 box_cell[0,0] = a_size
 box_cell[1,0] = b_size*np.cos(gamma)
 box_cell[1,1] = b_size*np.sin(gamma)
 box_cell[2,0] = c_size*np.cos(beta)
 box_cell[2,1] = c_size*( np.sin(beta)*( np.cos(alpha)-np.cos(beta)*np.cos(gamma) )/ \
                         (np.sin(beta)*np.sin(gamma)) )
 box_cell[2,2] = c_size*(np.sin(beta)*                                \
             np.sqrt(                                                \
                     ( np.sin(gamma)**2                              \
                      -np.cos(beta)**2 - np.cos(alpha)**2            \
                      +2.0*np.cos(alpha)*np.cos(beta)*np.cos(gamma) \
                     )                                               \
                     )/(np.sin(beta)*np.sin(gamma))                  \
                     )
 return box_cell.T;



#LAMMPS related ..............


def print_structures_lammps(nat_per_box,itype,lammps_cell,rcell,NameFile=None):
  xlo=0.0
  ylo=0.0
  zlo=0.0
  xhi = lammps_cell[0, 0]
  yhi = lammps_cell[1, 1]
  zhi = lammps_cell[2, 2]
  xy = lammps_cell[0, 1]
  xz = lammps_cell[0, 2]
  yz = lammps_cell[1, 2]

  if NameFile is None:
    fout = open("cube.lmp","w")
  else:
    fout=open(NameFile,"w")

  fout.write("# add your comment ... \n")
  fout.write(" \n")
  fout.write("%i atoms\n"%nat_per_box)
  ntypes=len(np.unique(itype))
  fout.write("%i atom types\n"%ntypes)
  fout.write(" \n")
  fout.write(" %22.16f  %22.16f   xlo xhi\n"%(xlo, xhi))
  fout.write(" %22.16f  %22.16f   ylo yhi\n"%(ylo, yhi))
  fout.write(" %22.16f  %22.16f   zlo zhi\n"%(zlo, zhi))
  fout.write(" %22.16f  %22.16f  %22.16f   xy xz yz\n"%(xy, xz,  yz))
  fout.write(" \n")
  fout.write("Atoms\n")
  fout.write(" \n")

  for i in range(nat_per_box):
    fout.write("%8i  %4i  %20.12f %20.12f %20.12f\n"%(i+1, itype[i], rcell[i,0],rcell[i,1],rcell[i,2]))
  fout.close()
  return



def write_lammps(nat_per_box, boxfin, itype, rcell, cell, outFile = None):


    #boxini and boxfin are in the form that we see the matrix in input:
    #                      a1x, a1y, a1z = a1 = A
    #     boxini    =      a2x, a2y, a2z = a2 = B
    #                      a3x, a3y, a3z = a3 = C
    #
    #for the calculus should be used box=boxini.T
    #                      a1x, a2x, a3x =
    #     box       =      a1y, a2y, a3y = a1 a2 a3 = A B C
    #                      a1z, a2z, a3z =
    #
    #box are in the format  box[0-2,:]=(A,B,C) first index cartesian projections,
    #                                          second index axis

    #x_cart (3,nat) =  box    (3 ,3) x x_cryst(3,nat)
    #x_cryst(3,nat) =  box^-1 (3 ,3) x x_cart (3,nat)

    tri_mat, beta_matrix = convert_cell_from_gin_to_lammps(boxfin)
    # beta_matrix = < e_i', e_j> e_j and e_i' define the reper of old and new axis, repectively.
    # tri_max is in box format
    if beta_matrix is None:
       beta_matrix=np.matlib.eye(3,dtype=float)

    globalv.beta_matrix=beta_matrix

    transd=np.zeros((3,3))

    #for i in range(3):
    #  for j in range(3):
    #    transd[i,j]=np.dot(tri_mat.T[:,i],boxini.T[:,j])
    #
    #these two are equivalent ...
    #
    #trans=np.matmul(tri_mat,boxini.T)

    rcell_new=np.matmul(beta_matrix,rcell.T).T
    tri_mat_inv=np.linalg.inv(tri_mat)
    cell_new = np.matmul(tri_mat_inv,rcell_new.T).T


    #this is just for testing ....
    #box_inv=np.linalg.inv(boxini.T)
    #cell_new = np.dot(box_inv,rcell.T).T
    #print_structures_ndm(nat_per_box,1,boxini,rcell,cell_new)

    print_structures_lammps(nat_per_box, itype,  tri_mat, rcell_new,outFile)
    return
