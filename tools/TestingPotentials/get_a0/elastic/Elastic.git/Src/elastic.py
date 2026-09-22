#!/usr/local/bin/python
#Purpose : Calculate C11,C12 and C44 of BCC crystal structure

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
from bulk_and_a0          import *
from qha_volume           import *
from read_input           import *
from tools                import *
from lammps               import *
from polar_deformation    import main_polar, test_vasp_polar, test_lammps_polar
from common               import *
#from pandas import *


usage = """

           Usage: python elastic.py typerun

           typerun -  'pre_build' (preparing files - for complex task qha_bulk)
                      'build'   (preparing files) ;
                      'run'    (running the prepared files) ;
                      'extract' ( extract the output) ;
                      'post_extract' (for qha_bulk, polar)
                      'all' - ALL of above :)
           all the other arguments being transmitted using input_elastic file
           (which should be present in the directory where the script is launched)


"""

#strains with constant volume deformation ...


# Main Program
#==========================================================================
# Default setting
#=====================================================================================
globalv.init()
unit_conversion1 = globalv.unit_conversion1

script=os.path.basename(os.path.normpath(sys.argv[0]))
globalv.script=script

if ( (len(sys.argv)> 2 )  ):
        print ("\n\tERROR: Wrong number of arguments !!!. Passing variables in argument is no longer possible ")
        print (usage)


inputf_ini="input_elastic.ini"
if not os.path.exists(inputf_ini):
   print(("elastic: The input file  %s doesn't exist. Create the file with the correct input."%(str(inputf_ini))))
   exit(0)
read_input_file(inputf_ini)
if len(sys.argv)==2:
      globalv.typerun=sys.argv[1]
mode=globalv.mode

structure=globalv.structure
typerun=globalv.typerun

nt_strain=globalv.nt_strain
a=globalv.a_guess
cal_elas=globalv.cal_elas
dynamics=globalv.dynamics

if cal_elas !='polar':
   if dynamics != 'min':
      print(('dynamics %s is not implemented for cal_elas %s'%(dynamics, cal_elas)))









based=os.getcwd()
globalv.based=based
globalv.root_dir=based
perc_strain_min=globalv.perc_strain_min
perc_strain_max=globalv.perc_strain_max

if (perc_strain_min > perc_strain_max):
  print ('prec_strain_min should be lower than perc_strain_max')
  print(('perc_strain_min = %s perc_strain_max = %s'%(perc_strain_min, perc_strain_max)))
  exit(0)

dt = (perc_strain_max - perc_strain_min)/float(nt_strain)
globalv.dt=dt

if (cal_elas=='qha_bulk') or (cal_elas=='qha_cxx' or cal_elas=='ana_bulk' or cal_elas=='ana_cxx'):
  perc_vol_min=globalv.perc_vol_min
  perc_vol_max=globalv.perc_vol_max
  nt_strain=globalv.nt_strain
  nt_temp=globalv.nt_temp
  temp_min=globalv.temp_min
  temp_max=globalv.temp_max
  a0=globalv.a_guess
  xp_vol , step_xp_vol=np.linspace(perc_vol_min,perc_vol_max,nt_strain,endpoint=True,retstep=True)
  xp_temp,step_xp_temp=np.linspace(temp_min,temp_max,nt_temp,endpoint=True,retstep=True)

  #xp_c11,step_xp_c11=np.linspace(perc_c11_min,perc_c11_max,nt_c11,endpoint=True,retstep=True)
  #xp_c44,step_xp_c44=np.linspace(perc_c44_min,perc_c44_max,nt_c44,endpoint=True,retstep=True)

  print(("V_min    V_max    nV  :"'{0:10.2f}'.format(perc_vol_min*100.0),  '{0:10.2f}'.format(perc_vol_max*100.0), '{0:5}'.format(nt_strain)))
  print(("T_min    T_max    nT  :"'{0:10.2f}'.format(temp_min),  '{0:10.2f}'.format(temp_max), '{0:5}'.format(nt_temp)))


E0= 1.0
V0 = 1.0
if (  (mode != "vasp")   and (mode != "ndm") and (mode != "phondy") and (mode != "lammps")  ):
  print("\n\tERROR: This argument should be vasp, pwscf, ndm, phondy, lammps")
  print(usage)
  exit(0)
if ((structure != "bcc") and (structure != "fcc") and (structure != "c15") and (structure !="a15") and (structure != "gin") ):
  print(("c",structure,"c"))
  print("\n\tERROR: The structure argument should be bcc, fcc, c15, a15,  or gin")
  print(usage)
  exit(0)
if ((typerun != "build") and (typerun!='test') and (typerun != "run") and (typerun != "all")  \
      and (typerun != "extract") and (typerun != "pre_build") and (typerun != 'post_extract') ):
  print("\n\tERROR: This argument should be pre_build, build, run, extract, post_extract or all")
  print( usage)
  exit(0)

if (typerun=='test'):
    itest_flag=0


    if (mode == 'vasp') and ( (cal_elas == 'qha_bulk') or (cal_elas == 'polar')) :
      itest_flag=1

    if (mode == 'lammps') and ((cal_elas == 'polar')) :
      itest_flag=1


    if (mode == 'ndm') and ((cal_elas == 'qha_bulk')) :
      itest_flag=1

    if itest_flag==0:
      print('test mode is not implemented for mode %s and cal_elas = %s'%(mode, cal_elas))
      print('--the only possiblities--')
      print('vasp      qha_bulk')
      print('vasp      polar')
      print('lammps    polar')
      exit(0)


# set yes or no we need for a pre building. When you do vacancy or other defects you certainly need pre_build.
# The structures
# should be minimized before the phonon spectrum should be computed.
# However when you perform elastic constants with the temperature no need for this pre_build mode.
# It can be replaced by the relaxation no or yes ?
pre_build_requested=False
if ((cal_elas=='qha_bulk' or cal_elas=='qha_cxx') and (mode=='vasp')):
  pre_build_requested=True
globalv.pre_build_requested=pre_build_requested

if typerun=='build' and pre_build_requested:
  print('WARNING: this cal_elas = %s require a pre_build mode. I hope you know what you doing.'%cal_elas)

if (typerun == 'build') or (typerun == 'all'):
   if  pre_build_requested==False:
     os.system("rm -rf bulk")
     os.system("rm -rf c11 c12 c44")
     os.system("rm -rf babulk")

if (typerun=='test'):

    if (mode=='vasp') and (cal_elas=='qha_bulk'):
       test_vasp_qha_bulk(based,'volume',xp_vol)

    if (mode=='vasp') and (cal_elas=='polar'):
       test_vasp_polar(based,'deformation')

    if (mode=='lammps') and (cal_elas=='polar'):
       test_lammps_polar(based,'deformation')

    if (mode=='ndm') and (cal_elas=='qha_bulk'):
       test_ndm_qha_bulk(based,'volume',xp_vol)




if ((mode == 'ndm') or (mode == 'vasp') or ( mode == 'pwscf') or (mode == 'phondy') or (mode=='lammps')):

 #debug print E0, V0, a0 11
 if (cal_elas=='polar'):
  if (structure !='gin'):
    print("cal_elas=polar is implemented only with structure=gin. Now structure is %s."%structure)
    exit(0)

  main_polar(typerun,based,a,structure)


 if ((cal_elas=='c11') or (cal_elas=='cxx') or (cal_elas=='c44')):
  E0,V0,a0 = get_E0V0(typerun,based,a)
 if ((cal_elas=='c11') or (cal_elas=='cxx')):
  k11 = get_cxx(based,"c11",a0,0.5,0.0,E0,V0)*unit_conversion1
  #old_version2 k12 = get_cxx(based,"c12",a0,0.5,0.0,E0,V0)*2.0/3.0*unit_conversion1
  k12=get_bulk_and_a0(based,"babulk",a0,0.5,1.0,E0)
  #old_version2 c12 = (k12 - k11)/3.0
  #old_version2 c11 = (k12 + 2.0*k11)/3.0

  c12 = (3.0*k12 - k11)/3.0
  c11 = (3.0*k12 + 2.0*k11)/3.0
 if ((cal_elas=='c44') or (cal_elas=='cxx')):
  #for old distortion: c44 = get_cxx(based,"c44",a0,2.0,1.0,E0,V0)/2*unit_conversion1
  c44 = get_cxx(based,"c44",a0,2.0,0.0,E0,V0)*2.0*unit_conversion1

 if (typerun == 'extract') or (typerun == 'all'):


  if ((cal_elas=='c11') or (cal_elas=='c44') or (cal_elas=='cxx')):
    fouto = open("responses.txt",'w')
    fout = open("cxx.txt",'w')
    if (cal_elas=='c11'):
      S = "C11\t%f \nC12\t%f \n"%(c11,c12)
    if (cal_elas=='cxx'):
      S = "C11\t%f \nC12\t%f \nC44\t%f \n"%(c11,c12,c44)
    if (cal_elas=='c44'):
      S = "C44\t%f\t \n"%(c44)
    print(S)
    fout.write(S)
    fout.close()

    fouto.write(S)
    fouto.close()


#elif ((mode == 'bvasp') or (mode =='bndm') or (mode=='bpwscf') or (mode =='bphondy')) or (mode=='blammps'):

if cal_elas=='bulk' :
 #print("debug1.......:", mode, cal_elas)
 E0,V0,a0 = get_E0V0(typerun,based,a)
 k12=get_bulk_and_a0(based,"babulk",a0,0.5,1.0,E0)

if cal_elas=='qha_bulk':
 if (globalv.nt_bulk_run==1):
  E0,V0,a0 = get_E0V0(typerun,based,a0,'bulk')
 if structure=='gin':
   a0=1
 E0=1
 gen_volume_qha(based,'volume',a0,E0,xp_vol,step_xp_vol,xp_temp,step_xp_temp, def_type=1)


if cal_elas=='qha_cxx':
 #E0,V0,a0 = get_E0V0(typerun,based,a)
 a0=1
 E0=1

 if (len(globalv.list_of_deformations) > 3):
    print("Maximum number of deformations in qha_cxx mode is 3.  Now is:", globalv.list_of_deformations)
    exit(0)



 if not (typerun != 'extract') or (typerun != 'post_extract'):
   for def_type in globalv.list_of_deformations:
     print(def_type, globalv.list_of_deformations)
     # is the deformation for B (t,t,t,0,0,0), former k12
     if (def_type==1):
       k11 = gen_volume_qha(based,'volume',a0,E0,xp_vol,step_xp_vol,xp_temp,step_xp_temp,def_type=1)

     # is the deformation for c11 (t,-t,1.0/(1.0-t**2)-1.0,0,0,0)
     if (def_type==2):
       k12 = gen_volume_qha(based,'deformation002',a0,E0,xp_vol,step_xp_vol,xp_temp,step_xp_temp,def_type=2)

     # is the deformation for c44
     if (def_type==3):
       c44 = gen_volume_qha(based,'deformation003',a0,E0,xp_vol,step_xp_vol,xp_temp,step_xp_temp,def_type=3)
 else:
    k11 = gen_volume_qha(based,'volume',a0,E0,xp_vol,step_xp_vol,xp_temp,step_xp_temp,def_type=1)

     # is the deformation for c11 (t,-t,1.0/(1.0-t**2)-1.0,0,0,0)
    k12 = gen_volume_qha(based,'deformation002',a0,E0,xp_vol,step_xp_vol,xp_temp,step_xp_temp,def_type=2)

     # is the deformation for c44
    c44 = gen_volume_qha(based,'deformation003',a0,E0,xp_vol,step_xp_vol,xp_temp,step_xp_temp,dfef_type=3)

    c12 = (3.0*k12 - k11)/3.0
    c11 = (3.0*k12 + 2.0*k11)/3.0






 #this is the future gen_grid_temperature_volume_qha(based,'grid_VT',a0,E0,xp_vol,step_xp_vol,xp_temp,step_xp_temp)
