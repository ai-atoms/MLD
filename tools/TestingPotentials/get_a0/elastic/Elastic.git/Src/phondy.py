import os
import re
import sys
import math
import numpy as np


def setup_phondy(path_dir):

  #eampot="/home/marinica/Flow/INPUT/eamtab.potin_W_EAM3"
  #eampot="/home/marinica/NDM_AREA51/JP_NDM/INPUT/eamtab.potin_Cu_MISHIN"
  #eampot="/home/marinica/NDM_AREA51/JP_NDM/INPUT/eamtab.potin_My_PERFECT_lot"
  #eampot="/home/marinica/NDM_AREA51/JP_NDM/INPUT/eamtab.potin_DD06"
  #eampot="/home/marinica/NDM_AREA51/JP_NDM/INPUT/eamtab.potin_M08"
  #eampot="/home/marinica/NDM_AREA51/JP_NDM/INPUT/eamtab.potin_M3EX_25p"

  dirPHONDY=os.environ['HOME']+'/FitFormation/Elastic/Phondy'
  atomic=dirPHONDY + '/data/atomic_fe'
  miniin=dirPHONDY + '/data/data_PHONDY_V_NORELAX.in'
  inputb=dirPHONDY + '/data/input.bash'

  if not os.path.exists(atomic):
   print("atomic file %s doesn't exist. Put the correct path in setup_phondy"%str(atomic))
   exit(0)
  if not os.path.exists(miniin):
   print("mini.in file %s doesn't exist. Put the correct path in setup_phondy"%str(miniin))
   exit(0)
  if not os.path.exists(inputb):
   print("input.bash file %s doesn't exist. Put the correct path in setup_phondy"%str(inputb))
   exit(0)


  os.system("cd %s"%(path_dir))
  os.system("rm -rf data results")
  os.system("mkdir data")
  os.system("mv geom.data data/")
  os.system("mkdir results")
  os.system("ln -s %s  input.bash"%(inputb))
  os.chdir(path_dir+'/data')

  os.system("ln -s %s  atomic"%(atomic))
  os.system("ln -s %s  mini.data"%(miniin))

  os.chdir(path_dir)


def run_phondy(outPHONDY):
  exerun=os.environ['HOME']+'/Flow/CEA_POT_CYR_COS/src/LAST_EXE/mini_PR.x'
  if not os.path.exists(exerun):
   print("Exe file %s doesn't exist. Put the correct path in run_phondy"%str(exerun))
   exit(0)


  os.system(" %s < input.bash > %s "%(exerun, outPHONDY))

  return

def get_energy_from_phondy(outPHONDY):
 word=[]
 f=open(outPHONDY,'r')
 for line in f:
   if re.match("(.*) energie totale (.*)", line):
     word=line.split()
     #print word[1]
 f.close()
 return float(word[2]);

def get_volume_from_phondy(outPHONDY):
 word=[]
 f=open(outPHONDY,'r')
 for line in f:
   if re.match("(.*) volume cell  (.*)", line):
     word=line.split()
     #print word[1]
 f.close()
 return float(word[2]);


def get_natom_from_phondy(outPHONDY):
 word=[]
 f=open(outPHONDY,'r')
 for line in f:
   if re.match("(.*) energie totale (.*)", line):
     word=line.split()
     ene = float(word[2])
 f.close()
 f=open(outPHONDY,'r')
 for line in f:
   if re.match("(.*) energie par atome (.*)", line):
     word=line.split()
     eneat = float(word[3])
 nat=int(round(ene/eneat))  #the best way to have the nint in python
 return float(nat);



def print_structures_phondy(nat_per_box,nrepeat,boxfin,itype,rcell,cell,outFile=None):
  if np.amax(itype) > 1:
    print('WARNING: print_structure_phondy doesnt handle more than one type')
    exit(0)

  if outFile is None:
    fout = open("geom.data","w")
  else:
    fout = open(outFile,"w")

  fout.write("%i\n"%(nat_per_box))
  fout.write("1.0\n")
  fout.write(" %22.16f  %22.16f  %22.16f\n"%(boxfin[0,0], boxfin[0,1],boxfin[0,2]))
  fout.write(" %22.16f  %22.16f  %22.16f\n"%(boxfin[1,0], boxfin[1,1],boxfin[1,2]))
  fout.write(" %22.16f  %22.16f  %22.16f\n"%(boxfin[2,0], boxfin[2,1],boxfin[2,2]))

  for i in range(nat_per_box):
    fout.write("%15.7f %15.7f %15.7f\n"%(rcell[i,0],rcell[i,1],rcell[i,2]))

  fout.close()
  return 
