import os
import re
import sys
import math
import globalv
import setup_param_ndm
import numpy as np
from lammps import write_lammps, read_lammps_data
from numpy.linalg import norm

# setut_ndm    ... only for minimizations ...
# setup_ndm_ph ... active for phondy calcuulation. Called from qha_volume.py

def setup_ndm(based,indir,cname):

  setup_param_ndm.init()
  os.chdir(based+'/'+indir)
  eampot=setup_param_ndm.eampot
  dinfile=setup_param_ndm.dinfile
  mlfile=setup_param_ndm.mlfile
  computer=globalv.cluster

  os.system("ln -s %s  eamtab.potin"%(eampot))
  os.system("ln -s %s  bso4_snap1_params.pot"%(eampot))
  os.system("cp -p %s  cube.din"%(dinfile))
  os.system("cp -p %s  cube.ml"%(mlfile))
  os.system("cp -p %s  ."%(eampot))

  if computer=='gatsby':
    generate_jsub_ndm_gatsby(indir,cname)
  if computer=='occigen':
    generate_jsub_ndm_occigen(indir,cname)
  if computer=='irene':
    generate_jsub_ndm_irene(indir,cname)


def setup_ndm_ph(based,indir,cname):


  setup_param_ndm.init()
  dinfile_ph=setup_param_ndm.dinfile_ph
  exerun_ph=setup_param_ndm.EXENDM_PH
  computer=globalv.cluster
  nprocs=setup_param_ndm.ntasks

  # for HA
  flmodes=setup_param_ndm.flmodes
  fldos  =setup_param_ndm.fldos
  inpmab =setup_param_ndm.inpmab

  os.chdir(based+'/'+indir)
  os.system("ln -s %s cube.phondy"%(inpmab))
  os.system("ln -s %s cube.phondy.ldos"%(fldos))
  os.system("ln -s %s cube.phondy.lmodes"%(flmodes))


  if (globalv.cal_elas=='qha_bulk'):
    if(globalv.qha_minmode=='lammps'):
      os.system("ln -s %s min.in"%(globalv.mininput_ph))
      os.system("ln -s %s in.lammps"%(globalv.phinput_ph))
      os.system("ln -s %s pot.fs"%(globalv.minpot01_ph))
      os.system("ln -s %s pot2.fs"%(globalv.minpot02_ph))
      if (globalv.qha_minmode=='lammps'):
        nat_per_box, boxini, itype, cell, rcell= read_gin (dirSTR='./', ginFile='cube.gin')
        write_lammps(nat_per_box, boxini, itype, rcell, cell, outFile='cube.lmp')


def run_ndm(outNDM):
  script=globalv.script
  setup_param_ndm.init()
  computer=globalv.cluster

  if computer=='gatsby':
    run_ndm_gatsby(outNDM)

  if computer=='occigen':
    run_ndm_occigen(outNDM)

  if computer=='irene':
    run_ndm_irene(outNDM)

  if computer=='my_computer':
    run_ndm_my_computer(outNDM)

  return



def run_ndm_my_computer(outNDM):
  script=globalv.script
  cal_elas=globalv.cal_elas
  setup_param_ndm.init()
  exerun=setup_param_ndm.EXENDM
  ntasks=setup_param_ndm.ntasks

  #debug print 'r_my_computer', exerun
  if cal_elas !='qha_bulk':
    if (globalv.debug > 0 ):
      print("minimization using ndm 1 .............")

    if ntasks > 1:
      os.system("mpirun -np %s %s > %s "%(ntasks, exerun, outNDM))
    else:
      os.system(" %s > %s "%(exerun, outNDM))


  if cal_elas=='qha_bulk':
    dinfile_ph=setup_param_ndm.dinfile_ph
    exerun_ph=setup_param_ndm.EXENDM_PH
    relaxation=globalv.relaxation
    # here is ndm - min ; ndm - ph
    #exerun=setup_param_ndm.EXEMIN_PH
    exerun=globalv.minexe_ph
    if (globalv.qha_minmode=='ndm'):
      if (globalv.debug > 0 ):
        print("debug: minimization using ndm 01 .............")
      dinfile_cg=setup_param_ndm.dinfile_cg
      sauve2gin=setup_param_ndm.SAUVE2GIN
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

    if (globalv.qha_minmode=='lammps'):
      if relaxation=='no':
        print (ntasks)
        print(exerun)
        os.system('mpirun -np %s %s < min.in > out.run'%(ntasks, exerun))
      if relaxation=='yes':
        if (globalv.debug > 0 ):
          print("debug: minimization using lammps ntasks exe 01 .............", ntasks, exerun)
          print(globalv.structure)
        if globalv.structure=='gin':
          nat_per_box, boxini, itype, cell, rcell = read_gin (dirSTR='./', ginFile='cube.gin')
          write_lammps(nat_per_box, boxini, itype, rcell, cell, outFile='cube.lmp')    
        os.system('mpirun -np %s %s < min.in > out.run'%(ntasks, exerun))
        print("debug: minimization using lammps 02 .............")
        nat_per_box, boxini, itype, cell, rcell = read_lammps_data(dirSTR='./',name_of_input_file='lammps.data')
        print("debug: minimization using lammps 03 .............")
        os.system("mv cube.lmp cube.lmp_ini")
        print("debug: minimization using lammps 04 .............") 
        os.system("mv cube.gin cube.gin_ini")
        nrepeat=1
        print_structures_ndm(nat_per_box,nrepeat,boxini,itype,rcell,cell)
        write_lammps(nat_per_box, boxini, itype, rcell, cell, outFile='cube.lmp')

    if globalv.qha_phmode=='phondy-ndm' or globalv.qha_phmode=='phondy-lammps':
      if (globalv.debug > 0 ):
        print("debug: phonons using phondy-ndm 02 .............", ntasks, exerun_ph)
      os.system(" cp -p %s cube.din"%(setup_param_ndm.dinfile_ph))
      os.system(" mpirun -np %s  %s > %s "%(ntasks, exerun_ph, 'out.run_ph'))

  return


def run_ndm_gatsby(outNDM):

#  if not os.path.exists(exerun):
#   print "Exe file %s doesn't exist. Put the correct path in setup_vasp"%str(exerun)

  os.system('pwd')
  os.system('qsub jsub')
  return

def run_ndm_occigen(outNDM):

#  if not os.path.exists(exerun):
#   print "Exe file %s doesn't exist. Put the correct path in setup_vasp"%str(exerun)

  os.system('pwd')
  os.system('sbatch  jsub')
  return


def run_ndm_irene(outNDM):

#  if not os.path.exists(exerun):
#   print "Exe file %s doesn't exist. Put the correct path in setup_vasp"%str(exerun)

  os.system('pwd')
  os.system('ccc_msub  jsub')
  return


def generate_jsub_ndm_occigen(indir,cname):
#this jsub for gatsby my friends ...
  script=globalv.script
  cal_elas=globalv.cal_elas
  setup_param_ndm.init()
  exerun=setup_param_ndm.EXENDM
  nodes=globalv.nodes
  nprocs=globalv.ntasks
  type_node=globalv.type_node

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
     print(("nodes:  %s"%nodes))
     print(("nprocs:  %s"%nprocs))
     exit(0)

  f=open('jsub','w')
#!/bin/sh
#SBATCH -J Fe_1C
##SBATCH --nodes=29                # debug:  number of nodes
##SBATCH --constraint=BDW28          # debug:  type of nodes
#SBATCH --constraint=[HSW24|BDW28]
##SBATCH --ntasks-per-node=28      # debug:  tasks /  nodes
#SBATCH --ntasks=720               # the number of MPI Tasks
#SBATCH --exclusive
#SBATCH --threads-per-core=1      #
#SBATCH --time  23:20:55         # Time limit: HH:MM:SS format
##SBATCH -p highmem              #  By default 2000 MB are reserved. More uncomment
#SBATCH --output  ofile.out   # stdout pf the job
#SBATCH --error   efile.err   # stderr of the job




  f.write('#!/bin/sh\n')
  f.write('#SBATCH -J %s%s\n'%(indir,cname))
  f.write('#SBATCH --nodes=%s               # debug:  number of nodes \n'%nodes)
  f.write('#SBATCH --constraint=[HSW24|BDW28]\n')
  f.write('#SBATCH --ntasks=%s              # the number of MPI Tasks \n'%nprocs)
  f.write('#SBATCH --exclusive              #  \n')
  f.write('#SBATCH --threads-per-core=1     #  \n')
  f.write('#SBATCH --time  00:20:00         # Time limit: HH:MM:SS format\n')
  f.write('##SBATCH -p highmem              #  By default 2000 MB are reserved. More uncomment\n')
  f.write('#SBATCH --output  ofile.out   # stdout pf the job\n')
  f.write('#SBATCH --error   efile.err   # stderr of the job\n')
  f.write('\n')
  f.write('\n')
  f.write('#module purge\n')
  f.write('module load intelmpi/2017.0.098\n')
  f.write('module load intel\n')
  f.write('module load  python\n')
  f.write('##---------------------\n')
  f.write('pwd\n')
  f.write('## SLURM_SUBMIT_DIR = the directory from which sbatch was invoked \n')
  f.write('## SLURM_JOB_ID  = job ID of the executing job \n')
  f.write('\n')
  f.write('module list\n')
  f.write('\n')
  f.write('suff=%s\n'%(indir))
  f.write('export nprocs=${SLURM_NTASKS}\n')
  f.write('echo ${nprocs} \n')
  f.write('\n')
  f.write('SUBMISSION=${SLURM_SUBMIT_DIR}\n')
  f.write('export OMP_NUM_THREADS=1\n')
  f.write('\n')
  f.write('\n')
  f.write('pwd\n')
  f.write('\n')
  f.write('\n')
  f.write('outputfile=${SUBMISSION}/w${suff}.out\n')
  f.write('ulimit -s unlimited\n')
  f.write('#export MPI_GROUP_MAX=1024\n')
  f.write('\n')
  if script!='qha_bulk':
    f.write('mpirun -np %s  %s  > out.run\n'%(nprocs, exerun))


  if cal_elas=='qha_bulk':
    exerun_ph=setup_param_ndm.EXENDM_PH
    exerun=setup_param_ndm.EXENDM
    sauve2gin=setup_param_ndm.SAUVE2GIN
    dinfile_cg=setup_param_ndm.dinfile_cg
    relaxation=globalv.relaxation

    f.write('export EXE=%s\n'%exerun_ph)
    f.write(" cp -p %s cube.din\n"%(dinfile_cg))
    f.write("mpirun -np %s %s > out.run_cg\n"%(nprocs, exerun))
    if relaxation=='no':
     f.write("mv out.run_cg out.run\n")
    elif relaxation=='yes':
     dinfile_cg2=setup_param_ndm.dinfile_cg2
     dinfile_tr=setup_param_ndm.dinfile_tr
     ELASTIC=setup_param_ndm.ELASTIC

     f.write("python %s/Src/get_imm.py cube.din > imm.tmp \n"%ELASTIC)
     f.write("%s cube.cout cube.gin < imm.tmp >/dev/null \n"%sauve2gin)

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


def generate_jsub_ndm_gatsby(indir,cname):
#this jsub for gatsby my friends ...

  script=globalv.script
  cal_elas=globalv.cal_elas
  setup_param_ndm.init()
  exerun=setup_param_ndm.EXENDM
  nodes=setup_param_ndm.nodes
  nprocs=setup_param_ndm.ntasks

  f=open('jsub','w')
  f.write('#!/bin/sh\n')
  if (nprocs==1):
    f.write('#PBS -q prod\n')
  if (nprocs > 1):
    f.write('#PBS -q prod_para\n')
  f.write('#PBS -l select=%s:ncpus=%s\n'%(nodes,nprocs))
  f.write('#PBS -l walltime=10:00:00\n')
  f.write('#PBS -N %s%s\n'%(indir,cname))
  f.write('#PBS -e error\n')
  f.write('#PBS -o output\n')
  f.write('\n')
  f.write('cat  $PBS_NODEFILE\n')
  f.write('cd $PBS_O_WORKDIR\n')
  f.write('pwd\n')
  f.write('\n')
  f.write('# ./create_gin_from_inp.sh > out.create\n')
  f.write('\n')
  f.write('cd $PBS_O_WORKDIR\n')
  f.write('\n')
  f.write('\n')
  f.write('#phonons calculation ...\n')
  if cal_elas !='qha_bulk':
    f.write('%s  > out.run\n'%(exerun))
    f.write('sleep 5\n')

  if  cal_elas=='qha_bulk':
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



def get_natom_from_ndm(outNDM):
 f=open(outNDM,'r')
 word=[]
 for line in f:
   if re.match(' nombre d atomes ', line):
     word=line.split()
     #print word[4]
 f.close()

 if len(word)==0:
   print(('get_natom_from_ndm: There is no nat output in the directory %s'%os.getcwd()))
   exit(0)

 return float(word[4]);

def get_energy_from_ndm(outNDM):
 f=open(outNDM,'r')
 word=[]
 for line in f:
   if re.match('(.*) energie (.*)', line):
     word=line.split()
     #print word
 f.close()

 if len(word)==0:
   print(('get_energy_from_ndm: There is no energie output in the directory %s'%os.getcwd()))
   exit(0)
 return float(word[-1]);

#def get_volume_from_ndm(outNDM):
# f=open(outNDM,'r')
# for line in f:
#   if re.match("(.*)volume=(.*)" , line):
#     word=line.split()
#     volume=str(word[3]).replace("D","E",1)
# f.close()
# return float(volume);


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

def dotProd(a,b):
# dor or scalar product in n dimension
 dimension = len(a)
 c = 0.0
 for i in range(dimension):
  c  += a[i]*b[i]
 return c



def get_sigma_from_ndm(outNDM):
  f=open(outNDM,'r')
  word=[]
  for line in f:
     if re.match("(.*)1 sigma total(.*)", line):
       word=line.split()
       sigma_1=[float(word[7]), float(word[8]),float(word[9])]
     if re.match("(.*)2 sigma total(.*)", line):
       word=line.split()
       sigma_2=[float(word[7]), float(word[8]),float(word[9])]
     if re.match("(.*)3 sigma total(.*)", line):
       word=line.split()
       sigma_3=[float(word[7]), float(word[8]),float(word[9])]


  sigma=np.mat((sigma_1, sigma_2, sigma_3))
  return sigma






def get_imm_from_din(dinNDM):
 f=open(dinNDM,'r')
 for line in f:
   if 'imm' in line:
    word=line.replace("="," ").split()
 f.close()
 return word[1];


def read_gin (dirSTR=None, ginFile=None):
     if dirSTR is None:
       root_dir=globalv.root_dir
       dirSTR=root_dir+'/Structure'

     if ginFile is None:
       sgin= dirSTR+'/'+ 'structure.gin'
     else:
       sgin= dirSTR+'/'+ ginFile

     if not os.path.exists(sgin):
        print(("Gin file %s doesn't exist. Put the correct path in setup_ndm and structure"%str(sgin)))
        exit(0)

     # should be readed in gin ...
     fgin=open(sgin,'r')
     lines_fgin=fgin.read().splitlines()
     icount=0
     for i in range(len(lines_fgin)):
      if not lines_fgin[i]=='':
        word=lines_fgin[i].split()
        if not (word[0]=='#'):
          icount=icount+1
          if icount==1:
            cell = [float(word[0]), float(word[1]), float(word[2]) ]
            cell_duplicate =np.array(cell)
          if icount==2:
             vbox_a= np.array([ float(word[0]),  float(word[1]),  float(word[2]) ])
          if icount==3:
             vbox_b= np.array([ float(word[0]),  float(word[1]),  float(word[2]) ])
          if icount==4:
             vbox_c= np.array([ float(word[0]),  float(word[1]),  float(word[2]) ])
             #boxini = np.mat (( vbox_a, vbox_b, vbox_c))
             #debug print vbox_a, np.shape(vbox_a)
             boxini = np.vstack ([ vbox_a, vbox_b, vbox_c])

             #                                              (A)
             #Please pay atttention of the fact that boxini=(B), has the same form as the input  file
             #                                              (C)
             #internal box should be box=(A,B,C)=boxini.T
             #print boxini.T[:,0]
          if icount==5:
             nat_per_box=int(word[0])
             cell=np.matrix(np.zeros((nat_per_box,3)))
             itype=np.array(np.zeros(nat_per_box))
          if (icount > 5):
             cell[icount-6,0]=float(word[0])
             cell[icount-6,1]=float(word[1])
             cell[icount-6,2]=float(word[2])
             itype[icount-6]=int(word[3])


     #rcell=np.zeros((nat_per_box,3))
     #for i in range(nat_per_box):
     #   for j in range(3):
     #    for k in range(3):
     #     rcell[i,j]=rcell[i,j]+cell[i,k]*boxini[k,j]

     box=boxini.T
     rcell = np.matmul(box, cell.T).T 

     return nat_per_box, boxini, itype, cell, rcell

def get_volume_from_ndm(outVASP):
 f=open('cube.gin','r')
 line=f.readlines()
 n1=int(line[0].split()[0])
 n2=int(line[0].split()[1])
 n3=int(line[0].split()[2])
 abox = []
 bbox = []
 cbox = []
 word=line[1].split()
 for i in range(3):
   abox.append(float(word[i]))
 word=line[2].split()
 for i in range(3):
   bbox.append(float(word[i]))
 word=line[3].split()
 for i in range(3):
   cbox.append(float(word[i]))

 vect = crossProd(bbox,cbox)
 volume=dotProd(abox,vect)*float(n1*n2*n3)
 f.close()
 return volume

def print_structures_ndm(nat_per_box,nrepeat,boxfin,itype,rcell,cell, outFile=None):
  if outFile==None:
    outFile='cube.gin'

  outname=outFile.split(".")[0]
  fout = open(outFile,"w")

  fout.write("%i %i %i\n"%(nrepeat,nrepeat,nrepeat))
  fout.write(" %22.16f  %22.16f  %22.16f\n"%(boxfin[0,0], boxfin[0,1],boxfin[0,2]))
  fout.write(" %22.16f  %22.16f  %22.16f\n"%(boxfin[1,0], boxfin[1,1],boxfin[1,2]))
  fout.write(" %22.16f  %22.16f  %22.16f\n"%(boxfin[2,0], boxfin[2,1],boxfin[2,2]))
  fout.write("%i\n"%(nat_per_box))


  for i in range(nat_per_box):
    fout.write("%15.7f %15.7f %15.7f  %4i\n"%(cell[i,0],cell[i,1],cell[i,2],itype[i]))
  fout.close()

  #writing name.in
  fout = open("name.in","w")
  fout.write(outname)
  fout.close()
  return



def generate_jsub_ndm_irene(indir,cname):
#this jsub for gatsby my friends ...
  script=globalv.script
  cal_elas=globalv.cal_elas
  setup_param_ndm.init()
  exerun=setup_param_ndm.EXENDM
  nodes=globalv.nodes
  nprocs=globalv.ntasks
  type_node=globalv.type_node

  if (globalv.type_node=='milan'):
     task_per_node=24
  elif (globalv.type_node=='rome'):
     task_per_node=24
  elif (globalv.type_node=='skylake'):
     task_per_node=24
  else:
     print("Unknown type of node. For irene it should be used only:")
     print("        milan or rome or skylake")
     exit(0)

  if (nprocs < 0) and (nodes > 0):
     nprocs=int(nodes)*task_per_node

  if (nodes < 0) and (nprocs > 0) :
     nodes = int(float(nprocs) / float(task_per_node)) + 1

  if (nodes < 0) and (nprocs < 0):
     print("At least one from nodes and nprocs should be greater than zero")
     print(("nodes:  %s"%nodes))
     print(("nprocs:  %s"%nprocs))
     exit(0)

  f=open('jsub','w')


  f.write('#!/bin/bash\n')
  f.write('#MSUB -r %s%s\n'%(cname[:3], indir[-4:]))
  f.write('#MSUB -n  %s\n'%nprocs)
  f.write('#MSUB -c %s\n'%globalv.cpu_per_task)

  f.write('#MSUB -q %s\n'%globalv.type_node)
  if (globalv.queue=='test'):
   f.write('#MSUB -Q test\n') 
  f.write('#MSUB -T %s\n'%globalv.runtime)
  f.write('#MSUB -o ofile.out\n')
  f.write('#MSUB -e ofile.err\n')
  f.write('#MSUB -A %s\n'%globalv.project)
  #f.write('#MSUB -A den\n')
  f.write('#MSUB -m scratch,work,store\n')


  f.write('##---------------------------\n')
  f.write('module purge\n')
  f.write('module load mpi/openmpi/4.0.5 scalapack/mkl/21.3.0\n')
  f.write('module load fortran/inteloneapi/21.4.0\n')
  f.write('module load gnu/12.2.0\n')
  f.write('module load cmake/3.22.2\n')
  f.write('##---------------------------\n')
  f.write('unset OMP_DISPLAY_ENV\n')
  f.write('unset OMP_NUM_THREADS\n')
  f.write('unset KMP_AFFINITY\n')
  f.write('unset OMP_PLACES\n')
  f.write('unset OMP_PROC_BIN\n')
  f.write('##---------------------------\n')
  
  
  f.write('pwd\n')
  f.write('suff=%s\n'%(indir))
  f.write('export nprocs=${BRIDGE_MSUB_NPROC}\n')
  f.write('echo ${nprocs} \n')
  f.write('\n')
  f.write('SUBMISSION=${BRIDGE_MSUB_PWD}\n')
  f.write('export OMP_NUM_THREADS=1\n')
  f.write('\n')
  f.write('pwd\n')
  f.write('\n')
  f.write('outputfile=${SUBMISSION}/w${suff}.out\n')
  f.write('ulimit -s unlimited\n')
  f.write('\n')
  if cal_elas!='qha_bulk':
    f.write('mpirun -np %s  %s  > out.run\n'%(nprocs, exerun))

  if cal_elas=='qha_bulk':
    exerun_ph=setup_param_ndm.EXENDM_PH
    exerun=setup_param_ndm.EXENDM
    sauve2gin=setup_param_ndm.SAUVE2GIN
    dinfile_cg=setup_param_ndm.dinfile_cg
    relaxation=globalv.relaxation
    if (globalv.qha_minmode=='ndm'):
      f.write('export EXE=%s\n'%exerun_ph)
      f.write(" cp -p %s cube.din\n"%(dinfile_cg))
      f.write("mpirun -np %s %s > out.run_cg\n"%(nprocs, exerun))
      if relaxation=='no':
        f.write("mv out.run_cg out.run\n")
      elif relaxation=='yes':
        dinfile_cg2=setup_param_ndm.dinfile_cg2
        dinfile_tr=setup_param_ndm.dinfile_tr
        ELASTIC=setup_param_ndm.ELASTIC
  
        f.write("python %s/Src/get_imm.py cube.din > imm.tmp \n"%ELASTIC)
        f.write("%s cube.cout cube.gin < imm.tmp >/dev/null \n"%sauve2gin)
  
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
        
    if(globalv.qha_minmode=='lammps'):
      if relaxation=='no':
        print (ntasks)
        print(exerun)
        f.write('mpirun -np %s %s < min.in > out.run\n'%(nprocs, exerun))
      if relaxation=='yes':
        if (globalv.debug > 0 ):
          f.write("#debug: minimization using plammps ntasks exe 01 .............: %s %s \n"%(nprocs,globalv.minexe_ph))
        if globalv.structure=='gin':
          nat_per_box, boxini, itype, cell, rcell = read_gin (dirSTR='./', ginFile='cube.gin')
          write_lammps(nat_per_box, boxini, itype, rcell, cell, outFile='cube.lmp')    

        f.write('mpirun -np %s %s < min.in > out.run\n'%(nprocs, globalv.minexe_ph))
        f.write("mv cube.lmp cube.lmp_ini\n")
        f.write("mv cube.gin cube.gin_ini\n")

        f.write("#debug: minimization using plammps 02 .............\n")
        f.write('rm -f cube.gin\n')
        f.write('python3 %s lammps.data cube.gin\n'%(globalv.lammps_to_gin))
        f.write('python3 %s cube.gin cube.lmp\n'%(globalv.gin_to_lammps))
        #nat_per_box, boxini, itype, cell, rcell = read_lammps_data(dirSTR='./',name_of_input_file='lammps.data')
        
        f.write(" \n")
        f.write("#debug: minimization using plammps 03 .............\n")
        nrepeat=1

    if globalv.qha_phmode=='phondy-ndm' or globalv.qha_phmode=='phondy-lammps':
      f.write(" cp -p %s cube.din\n"%(setup_param_ndm.dinfile_ph))
      f.write(" ccc_mprun  %s > %s "%(exerun_ph, 'out.run_ph\n'))

          

  f.close()
  return

