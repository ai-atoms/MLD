import os
import re
import sys
import math
import numpy as np
import globalv
from globalv import *
from ndm                  import *
from phondy               import *
from vasp                 import *
from pwscf                import *
from structures           import *




def get_E0V0 (typerun,based,a0,inpDir=None,ginFile=None,sigma=None):
 cal_elas=globalv.cal_elas
 mode=globalv.mode
 dynamics=globalv.dynamics
 os.chdir(based)

 if inpDir is None:
    dirrun="bulk"
 else: 
    dirrun=inpDir

 if (typerun == 'build') or ( typerun == 'all' ):
   if globalv.abinitio==0:
     os.system("rm -rf %s"%dirrun)
     os.mkdir(dirrun)
   os.chdir(dirrun)
 #debug print 'E0', mode, based
 if (typerun == 'build') or ( typerun == 'all' ):


# building structures files
  # it neeeds: gen_structure, strain_tensor, inpDir, mode, dynamics, ginfile, a0  
  strain_tensor=[0,0,0,0,0,0]
  write_structure(a0,a0,a0, strain_tensor,inpDir,ginFile)

# setup modes
  if ( mode == 'ndm') :
    setup_ndm(based,dirrun,dirrun)
    if cal_elas=='qha_bulk':
      setup_ndm_ph(based,dirrun,dirrun)
  if ( mode == 'phondy') :
    setup_phondy(based+'/'+dirrun)
  if ( mode == 'vasp') :
    setup_vasp(based,dirrun,dirrun)
    if (cal_elas=='qha_bulk'): 
       setup_vasp_ph(based,dirrun,dirrun)
  if ( mode == 'pwscf') :
    setup_pwscf(based,dirrun,dirrun)
  if ( mode == 'lammps') :
    setup_lammps(based,dirrun,dirrun)

 if (typerun == 'pre_build') or ( typerun == 'all' ):
    predirrun='pre'+dirrun 
    if glovalv.abinitio==0:
      os.system("rm -rf %s"%predirrun)
      os.mkdir(predirrun)
    os.chdir(predirrun)
    gen_structure([a0,a0,a0],strain_tensor=[0,0,0,0,0,0],ginFile=ginFile)
    if ( mode == 'vasp'):
       setup_vasp(based,predirrun,predirrun)
       run_vasp("out.run")

 if (typerun == 'run') or ( typerun == 'all' ):
   os.chdir(based + '/'+dirrun)
   run_md(mode,based,dirrun,Exceptions=None)

 if (typerun == 'build') or (typerun == 'run') or (typerun=='pre_build'):
  E0=-1
  V0=1
  a0=a0
  sigma0=1
  nat0=1
 if (typerun == 'extract') or ( typerun == 'all') or (typerun == 'post_extract') :
  os.chdir(based+'/'+dirrun)
  if (mode == 'ndm') :
   E0 = get_energy_from_ndm("out.run")
   V0 = get_volume_from_ndm("out.run")
   if sigma:
     sigma0=get_sigma_from_ndm("out.run")
     nat0=get_natom_from_ndm("out.run")

  if (mode == 'lammps') :
   if (dynamics=='neb' and inpDir[1:5] != 'bulk'):
     nat0, E0, V0, sigma0 = get_all_from_neb_lammps("out.run")
   else:
     E0 = get_energy_from_lammps("out.run",Patched=True)
     V0 = get_volume_from_lammps("out.run",Patched=True)
     if sigma:
       sigma0=get_sigma_from_lammps("out.run",Patched=True)
       nat0=get_natom_from_lammps("out.run",Patched=True)


  if (mode == 'phondy') :
   E0 = get_energy_from_phondy("out.run")
   V0 = get_volume_from_phondy("out.run")
   if sigma:
     print("not yet implemented in get_E0V0 for phondy")
     exit(0)

  if (mode == 'vasp') :
   E0 = get_energy_from_vasp_outcar("OUTCAR")
   V0 = get_volume_from_vasp_outcar("OUTCAR")
   if sigma:
     sigma0=get_sigma_from_vasp("OUTCAR",Patched=True)
     nat0=get_natom_from_vasp_outcar("OUTCAR")


  if (mode == 'pwscf') or (mode == 'bpwscf'):
   E0 = get_energy_from_pwscf_outcar("out.run")
   V0 = get_volume_from_pwscf_outcar("out.run")
   if sigma:
     print("not yet implemented in get_E0V0for pwscf")
     exit(0)

 if sigma is None:
    return (E0, V0, a0)
 else: 
    return (E0, V0, a0, nat0, sigma0)

def out_file_names(mode,dynamics):

     if mode=='lammps': 
         outFile_min='cube.lmp'
         outFile_deb='cube_deb.lmp'
         outFile_fin='cube_fin.lmp'
     elif mode=='vasp':
         outFile_min='POSCAR'
         outFile_deb='POSCAR_deb'
         outFile_fin='POSCAR_fin'
     elif mode=='pwscf':
         outFile_min='cube.pwscf'
         outFile_deb='cube_deb.pwscf'
         outFile_fin='cube_fin.pwscf'
     elif mode=='ndm':
         outFile_min='cube.gin'
         outFile_deb='cube_deb.gin'
         outFile_fin='cube_fin.gin'
     elif mode=='phondy':
         outFile_min='geom.data'
         outFile_deb='geom_deb.data'
         outFile_fin='geom_fin.data'
     else:
         print('no such dynamics')
         exit(0)
     if dynamics=='min':
        return outFile_min
     if dynamics=='neb':
        return outFile_deb, outFile_fin


def   write_structure(ax,ay,az, strain_tensor,inpDir,ginFile):
#  uppper layer for embeding gen_structure.

  dynamics=globalv.dynamics
  mode=globalv.mode

  if dynamics=='min':
     outFile_min = out_file_names(mode, dynamics)
     gen_structure([ax,ay,az],strain_tensor,ginFile=ginFile,outFile=outFile_min)

  if dynamics=='neb':
     if inpDir[1:5]=='bulk':
        outFile_min = out_file_names(mode, 'min')
        gen_structure([ax,ay,az],strain_tensor,ginFile=ginFile,outFile=outFile_min)
     else:
        outFile_deb, outFile_fin = out_file_names(mode, dynamics)
        gen_structure([ax,ay,az],strain_tensor,ginFile='structure_dfct_deb.gin',outFile=outFile_deb)
        nat_per_box, boxini, itype, rcell, cell = \
          gen_structure([ax,ay,az],strain_tensor,ginFile='structure_dfct_fin.gin',outFile=outFile_fin, generic=True)
        fn=open('cube_fin.disp','w')
        fn.write('%i\n'%nat_per_box)
        for i in range(nat_per_box):
            fn.write(' %i %.15f %.15f %.15f\n'%(i+1, rcell[i,0],rcell[i,1],rcell[i,2]))
        fn.close()
  return;



def   setup_calculation (dirBase, dirDeformation, dirLocal):

      mode=globalv.mode

      indir_b='ebulk_'+dirLocal

      os.chdir(dirBase)
      os.chdir(dirDeformation)
      indir=dirLocal
      tempd=dirBase+'/'+dirDeformation
      cname=dirDeformation

      if (mode == 'ndm') :
         setup_ndm(tempd,indir,'babulk')
      if (mode == 'lammps') :
       
         setup_lammps(tempd,indir,'babulk')
         if globalv.bulk_deform=='yes':
           setup_lammps(tempd,indir_b,'babulk')

      if (mode == 'vasp') :
         setup_vasp(tempd,indir,cname)
         if globalv.bulk_deform=='yes':
           setup_vasp(tempd,indir_b,cname)

      if (mode == 'phondy') :
         setup_phondy(tempd+'/'+indir)
      if (mode == 'pwscf') :
         setup_pwscf(tempd,indir,'babulk')



