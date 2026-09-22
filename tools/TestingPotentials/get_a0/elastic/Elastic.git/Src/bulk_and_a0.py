import globalv
import os
from structures import gen_structure
from ndm       import *
from lammps    import *
from phondy    import *
from vasp      import *
from pwscf     import *
from my_eosfit import *
from globalv   import run_md


def get_bulk_and_a0(based,cname,a0,a,b,c):
  mode=globalv.mode
  typerun=globalv.typerun
  perc_strain_min=globalv.perc_strain_min
  nt_strain=globalv.nt_strain
  dt=globalv.dt
#dmcm atoms_per_unit=globalv.atoms_per_unit
#Calculate elastic constants
  energy = []
  volume = []
  print("... %s"%cname)
  os.system("echo %s > breport.%s"%(cname,mode))
  os.chdir(based)
  if (typerun == 'build'  or ( typerun == 'all') ):
   os.mkdir("%s"%cname)
  nvalue=int(10000)
  for i in range(int(nt_strain)+1):
    t = perc_strain_min+i*dt
    #Generate lattice according to strain tensor
    indir=str(nvalue+i)
    ax,ay,az = a0,a0,a0


    #One by one version ....
    tempd=based+'/'+cname
    #debug print "tmpd  %s"%tempd
    os.chdir(tempd)
    if (typerun == 'build') or (typerun == 'all' ):
      os.mkdir("%s"%str(indir))
      os.chdir("%s"%str(indir))
      gen_structure([ax,ay,az],strain_tensor=[t,t,t,0,0,0])


      if (mode == 'ndm') :
         setup_ndm(tempd,indir,'babulk')
      if (mode == 'lammps') or (mode == 'blammps'):
         setup_lammps(tempd,indir,'babulk')
      if (mode == 'vasp') or (mode == 'bvasp'):
         #print("setup_vasp")
         setup_vasp(tempd,indir,'babulk')
      if (mode == 'pwscf') or (mode == 'bpwscf'):
         setup_pwscf(tempd,indir,'babulk')
      if (mode == 'phondy') or (mode == 'bphondy'):
         setup_phondy(tempd+'/'+indir)

    if (typerun == 'run' ) or (typerun == 'all') :
      os.chdir(tempd+'/'+indir)
      run_md(mode,tempd,indir,Exceptions=None)


 #   if mode  == 'run':
 #     os.system("%s > out.ndm 2>&1"%exendm)
    if (typerun == 'extract') or (typerun == 'all') :
      #debug print tempd, indir
      os.chdir(tempd+'/'+indir)
      if (mode=='ndm') :
       nat=get_natom_from_ndm("out.run")
       energy.append(get_energy_from_ndm("out.run")/float(nat))
       volume.append(get_volume_from_ndm("out.run")/float(nat))
      if (mode =='blammps' or mode=='lammps') :
       nat=get_natom_from_lammps("out.run",Patched=True)
       energy.append(get_energy_from_lammps("out.run",Patched=True)/float(nat))
       volume.append(get_volume_from_lammps("out.run",Patched=True)/float(nat))
      if (mode == 'bphondy' or mode=='phondy' ):
       nat=get_natom_from_phondy("out.run")
       energy.append(get_energy_from_phondy("out.run")/float(nat))
       volume.append(get_volume_from_phondy("out.run")/float(nat))
      if (mode == 'bvasp' or mode=='vasp' ):
       nat=get_natom_from_vasp_outcar("OUTCAR")
       #debug print("1111", nat)
       #debug print("2222", get_volume_from_vasp_outcar("OUTCAR"))
       #debug os.system('pwd') #debug
       energy.append(get_energy_from_vasp_outcar("OUTCAR")/float(nat))
       volume.append(get_volume_from_vasp_outcar("OUTCAR")/float(nat))
       #print(volume)
      if (mode == 'bpwscf' or mode=='pwscf' ):
       nat=get_natom_from_pwscf_outcar("out.run")
       energy.append(get_energy_from_pwscf_outcar("out.run")/float(nat))
       volume.append(get_volume_from_pwscf_outcar("out.run")/float(nat))
 #
  bulk_modules=0
  structure=globalv.structure
  if (typerun == 'extract') or (typerun == 'all') :
    os.chdir(based)
    fout = open("%s_volume_energ"%cname,"w")
    for i in range(int(nt_strain)+1):
      fout.write("%22.16f %22.16f\n"%(volume[i],energy[i]))
    fout.close()
    (volume0, bulk_modules, energy0)= fit_bulk_and_a0(cname)
    print("  a0(A)   V/at(A^3)    B       Ecoh")
    if structure=='bcc':
      nfac=2
    elif structure=='fcc':
     nfac=4
    elif structure=='c15':
     nfac=24
     #nfac=6
    elif structure=='a15':
     nfac=8
    elif structure=='gin':
     nfac=nat
    print((" %7.5f %8.5f %8.1f  %14.8f "%((volume0*float(nfac))**(1.0/3.0), volume0, bulk_modules, energy0)))

    fa = open("a0.txt",'w')
    fa.write("%7.5f\n"%((volume0*float(nfac))**(1.0/3.0)))
    fa.write("%s\n"%bulk_modules)
    fa.close()


  return bulk_modules;

def get_cxx(based,cname,a0,a,b,c,V0):
#Calculate elastic constants
  energy = []
  mode=globalv.mode
  typerun=globalv.typerun
  nt_strain=globalv.nt_strain
  perc_strain_min=globalv.perc_strain_min
  perc_strain_max=globalv.perc_strain_max
  dt=globalv.dt
  E0=c

  print("... %s"%cname)
  os.system("echo %s > report.%s"%(cname,mode))
  os.chdir(based)
  if (typerun == 'build'  or ( typerun == 'all') ):
   os.mkdir("%s"%cname)
  nvalue=int(10000)
  for i in range(int(nt_strain)+1):
    t = perc_strain_min+float(i)*dt
    #Generate lattice according to strain tensor
    indir=str(nvalue+i)
    ax,ay,az = a0,a0,a0


    #One by one version ....
    tempd=based+'/'+cname
    #debug print "tmpd  %s"%tempd
    os.chdir(tempd)
    if (typerun == 'build') or (typerun == 'all' ):
      os.mkdir("%s"%str(indir))
    os.chdir("%s"%str(indir))
    if (typerun == 'build') or (typerun == 'all') :
      if cname == 'c11':
         gen_structure([ax,ay,az],strain_tensor=[t,-t,1.0/(1.0-t**2)-1.0,0,0,0])
      elif cname == 'c12':
         gen_structure([ax,ay,az],strain_tensor=[t,t,t,0,0,0])
      elif cname == 'c44':
         #old distorsion gen_bcc_ndmgin([ax,ay,az],strain_tensor=[0,0,0,2*t,0,0])

         gen_structure([ax,ay,az],strain_tensor=[0,0,t**2/(4.0-t**2),0,0,t])
         #old distorsion gen_bcc_ndmgin([ax,ay,az],strain_tensor=[0,0,0,2*t,0,0])


      if mode == 'ndm':
         setup_ndm(tempd,indir,cname)
      if mode == 'lammps':
         setup_lammps(tempd,indir,cname)
      if mode == 'vasp':
         setup_vasp(tempd,indir,cname)
      if mode == 'pwscf':
         setup_pwscf(tempd,indir,cname)
      if mode == 'phondy':
         setup_phondy(tempd+'/'+indir)


    if (typerun == 'run' ) or (typerun == 'all') :
      os.chdir(tempd + '/'+indir)
      run_md(mode, tempd, indir,Exceptions=None)

    if (typerun == 'extract') or (typerun == 'all') :
     os.chdir(tempd + '/'+indir)
     if (mode == 'ndm'):
       energy.append(get_energy_from_ndm("out.run"))
     if (mode == 'lammps'):
       energy.append(get_energy_from_lammps("out.run",Patched=True))
     if (mode == 'phondy'):
       energy.append(get_energy_from_phondy("out.run"))
     if (mode == 'vasp'):
       energy.append(get_energy_from_vasp_outcar("OUTCAR"))
     if (mode == 'pwscf'):
       energy.append(get_energy_from_pwscf_outcar("out.run"))


  if (typerun == 'extract') or (typerun == 'all') :
    os.chdir(based)
    fout = open("%s_summary"%cname,"w")
    for i in range(int(nt_strain)+1):
      fout.write("%22.16f %22.16f\n"%(perc_strain_min+float(i)*dt,energy[i]))
    fout.close()

    #debug print a0, V0, E0
    gpfile = gen_gpfile_for_parabola_fit(a*V0,1.0,E0,cname)
    os.system("gnuplot %s > gnuplot.report 2>&1"%gpfile)
    a,b,c = get_fitted_param(a*V0,b,c,"fit.log")
    #print 'a b c', a, b ,c
    os.system("mv fit.log %s_fit.log"%cname)
    if os.path.exists('report.'+mode):
       os.system("mv report.%s %s_report.%s"%(mode,cname,mode))
  return a/V0
