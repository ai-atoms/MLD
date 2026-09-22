import sys
import os
import numpy as np
import globalv
from vasp   import *
from lammps import *

def delta(i,j):
    if i == j:
        return 1.0
    else:
        return 0.0

def delta4(i,j,k,l):
    if i == j and i == k and  i == l:
        return 1.0
    else:
        return 0.0 

def get_box_for_aneto(dirbased, dircalc):

  os.chdir(dirbased)
  os.chdir(dircalc)

  if globalv.mode=='vasp':
    volume, boxsim = get_volume_from_vasp_outcar(outVASP='OUTCAR', box=True)
  if globalv.mode=='lammps':
    nat_per_box, boxsim, itype, cell, rcell = read_lammps_data(dirSTR=dirbased+'/'+dircalc, name_of_input_file='cube.lmp')
  if globalv.mode=='ndm':
    nat_per_box, boxsim, itype, cell, rcell = read_gin (dirSTR=dirbased+'/'+dircalc, ginFile='cube.gin')

  volume = np.fabs(np.dot(np.cross(boxsim[0,:],boxsim[1,:]),boxsim[2,:]))
  os.chdir(dirbased)
  return volume, boxsim;


def write_input_aneto_cubic(dirbased, bulkdir, sigma_bulk, box_bulk, defectdir, sigma_defect, box_defect, C11,C12,C44):
  
  os.chdir(dirbased)
  os.chdir(defectdir)

  f1 = open("input_elast" , "w")
  f1.write("&input\n")
  f1.write("    CVoigt(1,1)= %23.16e \n" % C11)
  f1.write("    CVoigt(1,2)= %23.16e \n" % C12)
  f1.write("    CVoigt(1,3)= %23.16e \n" % C12)
  f1.write("    CVoigt(2,1)= %23.16e \n" % C12)
  f1.write("    CVoigt(2,2)= %23.16e \n" % C11)
  f1.write("    CVoigt(2,3)= %23.16e \n" % C12)
  f1.write("    CVoigt(3,1)= %23.16e \n" % C12)
  f1.write("    CVoigt(3,2)= %23.16e \n" % C12)
  f1.write("    CVoigt(3,3)= %23.16e \n" % C11)
  f1.write("    CVoigt(4,4)= %23.16e \n" % C44)
  f1.write("    CVoigt(5,5)= %23.16e \n" % C44)
  f1.write("    CVoigt(6,6)= %23.16e \n" % C44)
  f1.write("\n")
  stress=np.zeros((3,3))
  if (sigma_defect.shape==sigma_bulk.shape):
     if sigma_defect.shape==(3,3):
       stress = (sigma_defect - sigma_bulk)/10.0
     elif sigma_defect.shape==(1,6):
       st = (sigma_defect - sigma_bulk)/10.0
       stress[0,0]=st[0,0]
       stress[1,1]=st[0,1]
       stress[2,2]=st[0,2]
       stress[1,2]=st[0,3]
       stress[0,2]=st[0,4]
       stress[0,1]=st[0,5]
       stress[1,0]=stress[0,1]
       stress[2,0]=stress[0,2]
       stress[2,1]=stress[1,2]
     else:
       print('in  write_input_aneto_cubic no sigma shape known',sigma_defect.shape)
       exit(0)
  else:
     print('in  write_input_aneto_cubic no sigma shape problems',sigma_defect.shape, sigma_bulk.shape)
     exit(0)

  A1_ref = box_bulk[0,:]
  A2_ref = box_bulk[1,:]
  A3_ref = box_bulk[2,:]
  alat=1.0

  f1.write("    sigma_res(1,1)= %23.16e\n" % stress[0,0])
  f1.write("    sigma_res(1,2)= %23.16e\n" % stress[0,1])
  f1.write("    sigma_res(1,3)= %23.16e\n" % stress[0,2])
  f1.write("    sigma_res(2,1)= %23.16e\n" % stress[1,0])
  f1.write("    sigma_res(2,2)= %23.16e\n" % stress[1,1])
  f1.write("    sigma_res(2,3)= %23.16e\n" % stress[1,2])
  f1.write("    sigma_res(3,1)= %23.16e\n" % stress[2,0])
  f1.write("    sigma_res(3,2)= %23.16e\n" % stress[2,1])
  f1.write("    sigma_res(3,3)= %23.16e\n" % stress[2,2])
  f1.write("\n")
  f1.write("    A1_ref(1) = %23.16e \n" % A1_ref[0])
  f1.write("    A1_ref(2) = %23.16e \n" % A1_ref[1])
  f1.write("    A1_ref(3) = %23.16e \n" % A1_ref[2])
  f1.write("    A2_ref(1) = %23.16e \n" % A2_ref[0])
  f1.write("    A2_ref(2) = %23.16e \n" % A2_ref[1])
  f1.write("    A2_ref(3) = %23.16e \n" % A2_ref[2])
  f1.write("    A3_ref(1) = %23.16e \n" % A3_ref[0])
  f1.write("    A3_ref(2) = %23.16e \n" % A3_ref[1])
  f1.write("    A3_ref(3) = %23.16e \n" % A3_ref[2])
  f1.write("    alat = %23.16e\n" %alat)

  f1.write("&end\n")
  f1.close()
  os.chdir(dirbased)
  return



#Cijkl = np.zeros((3,3,3,3))
#for i in range(3):
#    for j in range(3):
#        for k in range(3):
#            for l in range(3):
#                Cijkl[i,j,k,l] = C12*delta(i,j)*delta(k,l) \
#                               + C44*(delta(i,k)*delta(j,l)+delta(i,l)*delta(j,k)) \
#                               + (C11-C12-2e0*C44)*delta4(i,j,k,l)
#Omega = np.zeros((3,3))
#for k in range(3):
#    for l in range(3):
#        Omega[k,l] = Omega1*unit_n[k]*unit_n[l] + Omega2/3e0*delta(k,l)
#
#Pij=np.zeros((3,3))
#for i in range(3):
#    for j in range(3):
#        for k in range(3):
#            for l in range(3):
#                Pij[i,j] += Cijkl[i,j,k,l]*Omega[k,l]


def get_correction_aneto(dirbased, dircalc, fileAneto):

  os.chdir(dirbased)
  os.chdir(dircalc)
  
  f2=open(fileAneto,'r')

  for line in f2:
    if "Interaction energy with periodic images:  Eint" in line:
       x=line.split()
       Ecorr=float(x[-2])/2.0

  f2.close()
  os.chdir(dirbased)

  return Ecorr



def detect_if_aneto_was_there (fileAneto):

  try:
     faneto=open(fileAneto)
     faneto.close()
     laneto=True
     if (os.stat(fileAneto)==0):
       laneto=False
  except IOError:
     laneto=False
  #debug print laneto, os.system('pwd')
  return laneto 


def run_aneto(dirbased, dircalc):

  os.chdir(dirbased)
  os.chdir(dircalc)
  os.system("%s > aneto.out"%globalv.aneto_exe)
  os.chdir(dirbased)

  return;

def pack_aneto (dirbased, bulkdir, sigma_bulk, defectdir,sigma_defect):

  if not(detect_if_aneto_was_there(fileAneto='aneto.out')):
    volume, box_bulk = get_box_for_aneto(dirbased=dirbased, dircalc=bulkdir)
    volume, box_dfct = get_box_for_aneto(dirbased=dirbased, dircalc=defectdir)
    write_input_aneto_cubic(dirbased=dirbased, \
                            bulkdir=bulkdir, sigma_bulk=sigma_bulk,box_bulk=box_bulk, \
                            defectdir=defectdir,sigma_defect=sigma_defect,box_defect=box_dfct,\
                            C11=globalv.valc11,C12=globalv.valc12, C44=globalv.valc44)

    run_aneto(dirbased=dirbased, dircalc=defectdir)
  ecorr=get_correction_aneto(dirbased=dirbased, dircalc=defectdir, fileAneto='aneto.out')
  return ecorr;
