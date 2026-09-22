import globalv
import os
import scipy as sp
from scipy.interpolate import splrep,splev, UnivariateSpline

from structures import gen_structure
from ndm       import *
from phondy    import *
from vasp      import *
from pwscf     import *
from lammps    import *
from my_eosfit import *
from tools     import *
from globalv   import run_md

def gen_volume_qha(based,cname,a0,E0,xp_vol,step_xp_vol,xp_temp,step_xp_temp,def_type=None):
  mode=globalv.mode
  pre_build_requested=globalv.pre_build_requested
  typerun=globalv.typerun
  debug=globalv.debug
  if mode=='vasp' or mode=='lammps':
    """
    #this is order to prevent large number of jobs ... for 1 ph calculations hundreds of force
    #calculations are engaged, being very likely to overcome the largest limit of submited jobs in the queue
    #(e.g. for occigen, the largest number is 150).
    """
    nt_strain_min=globalv.nt_strain_min
    nt_strain_max=globalv.nt_strain_max
  elif mode=='ndm':
    nt_strain_min=0
    nt_strain_max=len(xp_vol)
  else:
    print(('not yet implemented for this mode = %s'%mode))
    print('only ndm or vasp')
    exit(0)

  if (mode == 'ndm') :
     setup_param_ndm.init()
     limit_omega=globalv.limit_omega
     npointq=1
  if (mode == 'pwscf') :
     setup_param_pwscf.init()
     limit_omega=globalv.limit_omega
     prefix=setup_param_pwscf.prefix
  if (mode == 'vasp') :
     setup_param_vasp.init()
     limit_omega=globalv.limit_omega
     prefix=setup_param_vasp.prefix
     npointq=1


  if (mode == 'lammps') :
     setup_param_lammps.init()
     limit_omega=globalv.limit_omega
     prefix=setup_param_vasp.prefix
     npointq=1



#dmcm atoms_per_unit=globalv.atoms_per_unit
#Calculate elastic constants
  energy = []
  volume = []
  print(("... %s"%cname))
  os.system("echo %s > breport.%s"%(cname,mode))
  os.chdir(based)
  if (typerun == 'build'  or ( typerun == 'all') ):
    if pre_build_requested == False :
      os.system('rm -rf %s'%cname)
      os.mkdir("%s"%cname)
    else:
      if not os.path.exists(cname):
         print("Probably you need to run pre_build before.")
         print(("Now the typerun is %s and the directory %s is missing"%(str(typerun), str(cname))))
         #exit(0)

  if (typerun == 'pre_build'  or ( typerun == 'all') ):
      os.system('rm -rf %s'%cname)
      os.mkdir("%s"%cname)

  nvalue=int(10000)
  itest=0
  #debug for i in range(len(xp_vol)):
  #debug   print i, xp_vol[i]
  #debug exit(0)
  for i in range(len(xp_vol)):
    t = xp_vol[i]
    #Generate lattice according to strain tensor
    indir=str(nvalue+i)
    #One by one version ....
    tempd=based+'/'+cname
    #debug print "tmpd  %s"%tempd
    os.chdir(tempd)
    if (typerun == 'build') or (typerun == 'all' ):
      os.system('rm -rf %s'%str(indir))
      os.mkdir("%s"%str(indir))
      os.chdir("%s"%str(indir))

      strain=apply_strain_elastic_constants (def_type,t)
      gen_structure([a0,a0,a0],strain_tensor=strain)
      #gen_structure([a0,a0,a0],strain_tensor=[t,t,t,0,0,0])

      if (debug > 1):
        print("debug:  build set the input in qha_volume for the volume", i)
      if (mode == 'ndm') :
         setup_ndm(tempd,indir,'babulk')
         setup_ndm_ph(tempd,indir,'babulk')
      if (mode == 'vasp') :
         setup_vasp(tempd,indir,'babulk')
         setup_vasp_ph(tempd,indir,'babulk')

      if (mode == 'lammps') :
         setup_lammps(tempd,indir,'babulk')
         setup_lammps_ph(tempd,indir,'babulk')

      if (mode == 'phondy') :
         setup_phondy(tempd+'/'+indir)
      if (mode == 'pwscf') :
         setup_pwscf(tempd,indir,'babulk')

    if (typerun == 'pre_build') or (typerun == 'all' ):
      os.mkdir("%s"%'pre'+str(indir))
      os.chdir("%s"%'pre'+str(indir))

      strain=apply_strain_elastic_constants (def_type,t)
      gen_structure([a0,a0,a0],strain_tensor=strain)
      #gen_structure([a0,a0,a0],strain_tensor=[t,t,t,0,0,0])



      if (mode == 'vasp') :
         setup_vasp(tempd,'pre'+indir,cname)
         run_vasp("out.run")

      if (mode == 'lammps') :
         setup_lammps(tempd,'pre'+indir,cname)
         run_vasp("out.run")


    if (typerun == 'run' ) or (typerun == 'all') or (typerun=='pre_run'):
      os.chdir(tempd+'/'+indir)
      if (i >= nt_strain_min) and (i < nt_strain_max):
        run_md(mode,tempd,indir,Exceptions=['phondy','pwscf'])


    if (typerun == 'extract') or (typerun == 'all') or (typerun=='post_extract'):
      #debug print tempd, indir
      if (mode=='ndm') :
        os.chdir(tempd+'/'+indir)
        if (globalv.qha_minmode=='ndm'):
          nat=get_natom_from_ndm("out.run")
          energy.append(get_energy_from_ndm("out.run")/float(nat))
          volume.append(get_volume_from_ndm("out.run")/float(nat))
        if (globalv.qha_minmode=='lammps'):
          nat=get_natom_from_lammps("out.run", Patched=True)
          energy.append(get_energy_from_lammps("out.run", Patched=True)/float(nat))
          volume.append(get_volume_from_lammps("out.run", Patched=True)/float(nat))

        print((" qha_volume read %s "%indir))
        omega=read_eigenvalues_from_ndm('eigenvalues.dat')
        #print(" ... done\n")
        for ii in  range(3,omega.size):
          if math.fabs(omega[ii]) >= limit_omega:
            if (omega[ii]<0.0):
              itest=1
              print(('WARNING imaginary freq for vol omega', i, omega[ii]))

      if (mode=='pwscf') :
        os.chdir(tempd+'/'+indir)
        nat=get_natom_from_pwscf_outcar("out.run")
        energy.append(get_energy_from_pwscf_outcar("out.run")/float(nat))
        volume.append(get_volume_from_pwscf_outcar("out.run")/float(nat))
        npointq, omega=read_eigenvalues_from_pwscf('%s.matdyn.modes'%prefix)
        for ii in  range(3,omega.size):
          if math.fabs(omega[ii]) >= limit_omega:
            if (omega[ii]<0.0):
              itest=1
              print(('WARNING imaginary freq for vol omega', i, omega[ii]))


      if (mode=='phondy'):
        print('qha_volume, phondy: this mode not yet implemented')
        exit(0)
        nat=get_natom_from_phondy("out.run")
        energy.append(get_energy_from_phondy("out.run")/float(nat))
        volume.append(get_volume_from_phondy("out.run")/float(nat))

      if (mode=='vasp'):
        #going in pre_directories
        os.chdir(tempd+'/'+'pre'+indir)
        nat=get_natom_from_vasp_outcar("OUTCAR")
        energy.append(get_energy_from_vasp_outcar("OUTCAR")/float(nat))
        volume.append(get_volume_from_vasp_outcar("OUTCAR")/float(nat))


        # going in indir directories
        os.chdir(tempd+'/'+indir)
        print((" read %s "%indir))
        #goes in each mresults_xxx directory
        if (typerun=='extract'):
         for ind in range(len(globalv.list_of_disp)):
           disp=globalv.list_of_disp[ind]
           dir_mresults='mresults_'+str(int(disp))
           os.chdir(dir_mresults)
           with open('tmpstr/list') as flist:
             fdisp=flist.read().splitlines()
           fforces=open('FORCES','w')
           fforces.write('%s\n'%(len(fdisp)))
           fforces.close()

           for tmp in fdisp:
              tmp_dir=tmp.strip()
              #this part should be removed after correcting old calculs... from here
              #tt bohr=0.5291772
              #tt with open('tmpstr/%s.inp'%tmp_dir) as finp:
              #tt   elem=finp.read().splitlines()
              #tt if len(elem) > 2:
              #tt   print 'problems in reading tmpstr', len(elem)
              #tt   print elem[0]
              #tt ytemp=elem[0].split()
              #tt fout=open('FORCES','a')
              #tt fout.write("  %s  %s  %s  %s\n"%(ytemp[0], '{0:14.9f}'.format(float(ytemp[1])/bohr), \
              #tt                              '{0:14.9f}'.format(float(ytemp[2])/bohr), \
              #tt                              '{0:14.9f}'.format(float(ytemp[3])/bohr) ))
              #tt fout.close()
              #tt #...........to here

              #if the tt patch are activated the next line should be comment.
              os.system('cat tmpstr/%s.inp >> FORCES'%tmp_dir)
              os.chdir(tmp_dir)
              if not os.path.exists('OUTCAR'):
                  print("vasp  file OUTCAR  doesn't exist in directory  ")
                  print((tempd+'/'+indir+'/'+dir_mresults+'/'+tmp_dir))

              iout = test_if_vasp_terminated_correctly("OUTCAR")
              if iout == 0:
                  print('WARNING the displacements is not computed. The job should be relanced')
                  print((tempd+'/'+indir+'/'+dir_mresults+'/'+tmp_dir))

              #get forces from all the list:
              nionic = get_nionic_step_from_vasp_outcar("OUTCAR")
              print_last_forces_from_vasp("OUTCAR",tempd+'/'+indir+'/'+dir_mresults+'/FORCES',nat,nionic)
              os.chdir('../')

           os.system('cp %s  INPHON'%(setup_param_vasp.dinfile_ph_02))
           search_word_replace_line('INPHON', 'DISP', 'DISP='+str(disp) +'\n')
           os.system('%s > /dev/null'%setup_param_vasp.EXENDM_PH)
           omega = read_eigenvalues_from_vasp('aFREQ.all')

        elif (typerun=='post_extract'):
          omega = read_eigenvalues_from_vasp('mresults_200/aFREQ.all')

        for ii in  range(3,omega.size):
          if math.fabs(omega[ii]) >= limit_omega:
            if (omega[ii]<0.0):
              itest=1
              print(('WARNING imaginary freq for vol omega', i, omega[ii]))



      #we pack all the omega for various volumes in
      #omega_all[i_volume,j_frequency_of_the_mode]
      # which is the jth frequency for the  ith volume
      if (i==0):
        omega_all=omega
      if (i>0):
        omega_all=np.vstack([omega_all,omega])
      #print i, np.shape(omega_all)
  #everything was readed ...


  #here is the  extract ...
  if (typerun == 'extract') or  (typerun == 'all') or (typerun=='post_extract') :
    if itest==1:
       print('extract is stopped in qha_volume cause of negative frequencies')
       exit(0)
    os.chdir(based)
    fab=open("out_a0_Bxx.dat","w")
    fout=open('out.dat','w')
    fab.write("#   T         a0           V0            B^T        B^S       Cv         Cp    alpha_V          alpha_l         Free \n")
    print("#   T         a0           V0            B^T        B^S       Cv         Cp    alpha_V            alpha_l          Free  \n")

    #Method 1 at each temperature we get the equilibrium value
    norder=4
    vol_min=volume[0]
    vol_max=volume[-1]

    test_vol , step_test_vol=np.linspace(vol_min,vol_max,1000,endpoint=True,retstep=True)
    Bxx=[]
    Bxxtemp=[]
    axx=[]
    Cvx=[]
    vxx=[]
    Fxx=[]
    for it in range(len(xp_temp)):
      yp=[]
      syp=[]
      #print 'temperature', it, xp_temp[it],len(xp_vol)
      #debug ft=open('temp%itest'%it,'w')
      #debug print (' temp %s'%xp_temp[it])
      print("SHAPE", omega_all.shape)
      for iv in range(len(xp_vol)):
         if iv ==0:
           fo,so,dso, uo = free_energy_classical(omega[3:],limit_omega,xp_temp[it])   
         if iv != 0:
           fo,so,dso, uo = free_energy_classical(omega_all[iv,3:],limit_omega,xp_temp[it])
            
         #fo,so,dso, uo = free_energy_quantum  (omega_all[iv,3:],limit_omega,xp_temp[it])
         #fo,so,dso, uo = free_energy_classical_old_slow  (omega_all[iv,:],limit_omega,xp_temp[it])
         yp.append(energy[iv]*float(nat) + fo/float(npointq))
         #debug print('volume npointq  iv %i %i %14.6f %14.6f %14.6f'%(npointq, iv,energy[iv]*float(nat),fo/float(npointq), energy[iv]*float(nat)+fo/float(npointq)))
         #debug ft.write('%6.1f %14.6f   %14.6f  %14.6f  \n'%(iv,energy[iv]*float(nat),fo,fo+energy[iv]*float(nat)))
         syp.append(dso/(float(npointq)*float(nat)))
         
      

      if iv == 0: 
        my_temp = [10, 200, 400, 600]
        for temp in my_temp:
          fo,so,dso, uo = free_energy_classical(omega[3:],limit_omega, temp) 
          print("{:<10} {:<10} {:<10} {:<10} {:<10}".format(temp, fo, so, dso, uo))
      if iv == 0:  
        print("Only one volume. Nothig to fit. Nothing to extrapolate. Exit. ")
        exit(0)  
      #debug ft.close()
      # yp     - free energy
      # volume -
      evtoGpa=160.20506
      min_energy=float(min(energy)*float(nat))
      for i in range(len(yp)):
         yp[i] = yp[i] - min_energy
         #print(("%6.1f  %14.6f   %14.6f  %14.6f  %14.6f"%(xp_vol[i], volume[i], yp[i], energy[i]*float(nat), float(min_energy))))
      aval=np.polyfit(volume,yp,norder)
      s_aval=np.polyfit(volume,syp,norder)
      #print 'aval',len(aval), aval
      #first derivative
      new_aval=-np.polyder(aval,1)
      #second derivative
      new_aval2= np.polyder(aval,2)

      spl=UnivariateSpline(volume, yp, k=4, s=0)
      spld=spl.derivative(n=1)
      spld2=spl.derivative(n=2)
      sol_full=spl.derivative().roots()
      y_eval=spl(sol_full)
      #print('debug', sol_full, y_eval)
      a0val=sol_full[np.argmin(y_eval)]

      if (debug > 3):
        ft=open('int_%s.dat'%str(int(xp_temp[it])),'w')
        fs=open('val_%s.dat'%str(int(xp_temp[it])),'w')
        """
        yv=np.polyval(aval,test_vol)
        yvd=np.polyval(new_aval,test_vol)

        """
        yv=spl(test_vol)
        yvd=spld(test_vol)

        for ivv in range(len(test_vol)):
           ft.write("%17.5f %17.5f %17.5f\n"%(test_vol[ivv],yv[ivv], yvd[ivv]))
        for ivv in range(len(yp)):
           fs.write("%17.5f %17.5f\n"%(volume[ivv],yp[ivv],))
        ft.close()
        fs.close()

      # or you can do this trick ...
      ##for iloc in range(norder):
      ##   new_aval.append(-aval[iloc]*float(norder-iloc))
      #a0val,sol  = get_roots_of_the_polynom_1D(vol_min-0.1,vol_max+1.0,new_aval, xp_temp[it])

      if np.polyval(new_aval2,a0val) < 0 :
            print('BIG PROBLEM in solution)')
      # if len(sol) > 1:
      #     for itt in range(len(sol)):
      #        print 'deriv', np.polyval(new_aval2,sol[itt])


      #debug print 'debug qha', xp_temp[it], a0val,(a0val*2.0)**0.33333333
      """
      Free0val=np.polyval(aval,a0val)
      Btemp=a0val*np.polyval(new_aval2,a0val)/float(nat)
      Cv0val= np.polyval(s_aval,a0val)

      """
      Free0val=spl(a0val)
      Btemp=a0val*spld2(a0val)/float(nat)
      Cv0val= spl(a0val)

      B0val=Btemp*evtoGpa
      Cvx.append(Cv0val)
      Bxx.append(B0val)
      Bxxtemp.append(Btemp)
      structure=globalv.structure
      if structure=='bcc':
        nfac=2
      elif structure=='fcc':
        nfac=4
      elif structure=='c15':
       nfac=24
      elif structure=='c15':
       nfac=8
      elif structure=='gin':
       nfac=2

      axx.append((a0val*nfac)**(1.0/3.0))
      vxx.append(a0val)
      Fxx.append(Free0val)
    #uncomment that if you want polynom technique in order to have the derivatives of interpolated points.
    ## we will use spline interpolation instead: is more accurate for dense grid.
    #old pol new_order=5
    #old pol if len(xp_temp) > 10:
    #old pol    new_order=50
    #old pol if len(xp_temp) < 10:
    #old pol    new_order=len(xp_temp)/2
    #old pol apol=np.polyfit(xp_temp,vxx,new_order)
    #old pol dapol=np.polyder(apol,1)
    #volumic thermal expansion ...
    vol_slpd=sp.interpolate.splrep(xp_temp,vxx)
    vol_derv=sp.interpolate.splev(xp_temp,vol_slpd,1)
    #linear thermal expansion ...
    lin_slpd=sp.interpolate.splrep(xp_temp,axx)
    lin_derv=sp.interpolate.splev(xp_temp,lin_slpd,1)
    jmol2ev=1.036410*1.e-5
    for i in range(len(xp_temp)):

      #old pol alpha=np.polyval(dapol,xp_temp[i])/vxx[i]
      alpha_vol=vol_derv[i]/vxx[i]
      alpha_lin=lin_derv[i]/axx[i]
      # to have alpha_lin from alpha_v you do not neeed all this thing
      # alpha_v = 3 alpha_l - is just simple as that
      #(just from definitions and partial derivatives definitions)
      Cpx=Cvx[i]+xp_temp[i]*vxx[i]*alpha_vol**2*Bxx[i]/evtoGpa
      Bs = Bxx[i]*(1.0+ xp_temp[i]*vxx[i]*alpha_vol**2/Cvx[i])
      #debug print alpha_vol,xp_temp[i],vxx[i],Cvx[i]
      print(("%6.1f  %14.6f   %14.6f  %14.6f  %14.6f   %14.6f  %14.6f %14.6f %14.6f %14.6f"%(xp_temp[i], axx[i], vxx[i], Bxx[i], Bs, Cvx[i]/jmol2ev, Cpx/jmol2ev, 1.e+5*alpha_vol, 1.e+5*alpha_lin,float(Fxx[i])+float(min_energy))))
      fab.write("%6.1f  %14.6f   %14.6f  %14.6f  %14.6f   %14.6f  %14.6f %14.6f %14.6f %14.6f \n"%(xp_temp[i], axx[i], vxx[i], Bxx[i], Bs, Cvx[i]/jmol2ev, Cpx/jmol2ev, 1.e+5*alpha_vol, 1.e+5*alpha_lin, float(Fxx[i])+min_energy))
    fab.close()

  if ((typerun=='extract') or (typerun=='post_extract')) and (def_type==1):
    return Bxx, Bs;
  else:
    return;




def  test_ndm_qha_bulk(based,cname,xp_vol):
# based is the root of the script
# cname is the "volume", "deformation" etc
  mode=globalv.mode
  typerun=globalv.typerun
  setup_param_ndm.init()
  limit_omega=globalv.limit_omega


  #testing pre files
  nvalue=int(10000)
  for i in range(len(xp_vol)):
    indir=str(nvalue+i)
    #One by one version ....
    tempd=based+'/'+cname
    os.chdir(tempd)
    os.chdir("%s"%(str(indir)))
    if os.path.exists('eigenvalues.dat'):
      iout = test_if_ndm_terminated_correctly("eigenvalues.dat")
      if iout == 0:
        print(('eigenvalues unfinished: %s'%(tempd+'/'+cname+'/pre'+str(indir))))
        os.system('sbatch jsub')
    else:
        print(("eigenvalues missing     : %s"%(tempd+'/'+cname+'/pre'+str(indir))))



  return ;




def  test_vasp_qha_bulk(based,cname,xp_vol):
# based is the root of the script
# cname is the "volume", "deformation" etc
  mode=globalv.mode
  typerun=globalv.typerun
  setup_param_vasp.init()
  limit_omega=globalv.limit_omega
  prefix=setup_param_vasp.prefix


  #testing pre files
  nvalue=int(10000)
  for i in range(len(xp_vol)):
    indir=str(nvalue+i)
    #One by one version ....
    tempd=based+'/'+cname
    os.chdir(tempd)
    os.chdir("%s"%('pre'+str(indir)))

    if os.path.exists('OUTCAR'):
      iout = test_if_vasp_terminated_correctly("OUTCAR")
      if iout == 0:
        print(('OUTCAR unfinished: %s'%(tempd+'/'+cname+'/pre'+str(indir))))
    else:
        print(("OUTCAR missing     : %s"%(tempd+'/'+cname+'/pre'+str(indir))))

  #testing mresults files
  nvalue=int(10000)
  for i in range(len(xp_vol)):
    indir=str(nvalue+i)
    #One by one version ....
    tempd=based+'/'+cname
    os.chdir(tempd)
    #testing pre files
    os.chdir("%s"%(str(indir)))

    for ind in range(len(globalv.list_of_disp)):
      disp=globalv.list_of_disp[ind]
      dir_mresults='mresults_'+str(int(disp))
      os.chdir(dir_mresults)
      with open('tmpstr/list') as flist:
          fdisp=flist.read().splitlines()
      for tmp in fdisp:
         tmp_dir=tmp.strip()
         os.chdir(tmp_dir)
         if  os.path.exists('OUTCAR'):
           iout = test_if_vasp_terminated_correctly("OUTCAR")
           if iout == 0:
              print(('OUTCAR unfinished: %s'%(tempd+'/'+indir+'/'+dir_mresults+'/'+tmp_dir)))
              #run_vasp('out.run')
         else:
           print(("OUTCAR missing   : %s"%(tempd+'/'+indir+'/'+dir_mresults+'/'+tmp_dir)))
         os.chdir('../')


  return ;
def apply_strain_elastic_constants (def_type,t):

      if def_type == 1:
        strain_tensor=[t,t,t,0,0,0]
      elif def_type == 2:
         strain_tensor=[t,-t,1.0/(1.0-t**2)-1.0,0,0,0]
      elif def_type == 3:
         #old distorsion gen_bcc_ndmgin([ax,ay,az],strain_tensor=[0,0,0,2*t,0,0])
         strain_tensor=[0,0,t**2/(4.0-t**2),0,0,t]
         #old distorsion gen_bcc_ndmgin([ax,ay,az],strain_tensor=[0,0,0,2*t,0,0])
      else:
           print(("Deformation type not known. def_type=%s"%def_type))

      return strain_tensor;
