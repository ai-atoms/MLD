import matplotlib
#matplotlib.use('TkAgg')
matplotlib.use('Agg')
from matplotlib import pyplot as pl
import globalv
from globalv import run_md
import os
import scipy as sp
from scipy.interpolate import splrep,splev,interp1d,CubicSpline
from structures import gen_structure
from ndm       import *
from phondy    import *
from vasp      import *
from pwscf     import *
from my_eosfit import *
from tools     import *
from voigt     import *
from lammps    import *
from common   import get_E0V0, write_structure, setup_calculation
from aneto    import *
from polar_fit import *


def main_polar(typerun,based,a,structure):
  unit_conversion1=globalv.unit_conversion1
  bulk_deform=globalv.bulk_deform

  E0b,V0b,a0b,nat0b,sigma0b = get_E0V0(typerun,based,a,inpDir='0bulk',ginFile='structure_bulk.gin',sigma=True)
  E0d,V0d,a0d,nat0d,sigma0d = get_E0V0(typerun,based,a,inpDir='0dfct',ginFile='structure_dfct.gin',sigma=True)
  E0=E0b
  if (typerun=='extract') or (typerun=='all'):

    if (globalv.aneto=='yes') :
      ecorr0 = pack_aneto (dirbased=based, bulkdir='0bulk', sigma_bulk=sigma0b, defectdir='0dfct',sigma_defect=sigma0d)

    vsigma_b0=pack_sigma_in_vsigma(sigma0b)
    vsigma_d0=pack_sigma_in_vsigma(sigma0d)
    #debug printMatrixE(vsigma_d0)
    #debug printMatrixE(vsigma_b0)
    print("Polarizability tensor for the structure provided in gin file \n")

    Efor=formation_energy(nat0b,nat0d,E0b,E0d)
    print("P_ij^0: 0 order of dipole tensor in Voigt notation (eV):")
    Pij0=(vsigma_d0/10.0-vsigma_b0/10.0)*V0b/unit_conversion1
    printMatrixE(Pij0)
    print("Formation energy  (eV): ")
    print(Efor)


  nt_strain=globalv.nt_strain
  def_xp_strain=np.array([]).reshape(0,nt_strain)
  index_def=[]

  for i in globalv.list_of_deformations:

    perc_strain_min=globalv.perc_strain_min
    perc_strain_max=globalv.perc_strain_max
    if globalv.list_of_deformations_nostandard!="":
      for j in  range(len(globalv.index_list_of_deformations_nostandard)):
        if globalv.index_list_of_deformations_nostandard[j]==i:
          perc_strain_min=-globalv.value_list_of_deformations_nostandard[j]
          perc_strain_max= globalv.value_list_of_deformations_nostandard[j]
    xp_strain , step_xp_strain=np.linspace(perc_strain_min,perc_strain_max,nt_strain,endpoint=True,retstep=True)
    def_xp_strain=np.vstack([def_xp_strain,xp_strain])
#  for i in range(7):
#     print def_xp_strain[i,:]
#this should be change deformation after deformation ...

  sigma_tot=np.array([]).reshape(0,6)
  sigma_tot_b=np.array([]).reshape(0,6)
  strain_tot=np.array([]).reshape(0,6)
  ene_tot=np.array([]).reshape(0,1)
  ene_tot_corr=np.array([]).reshape(0,1)
  ene_tot_b=np.array([]).reshape(0,1)
  ene_int_tot=np.array([]).reshape(0,1)

  ii=0
  for def_type in globalv.list_of_deformations:
     ene_i, ene_i_b, ene_corr_i, sigma_i,sigma_i_b, strain_i= \
     get_polar_deformation(based,"deformation",def_type,a0d,E0,sigma0b, sigma0d,def_xp_strain[ii,:])
     ii+=1
     if (typerun=='extract') or (typerun=='all'):
       edipole=np.matmul(-Pij0,strain_i.T).T
       #debug print strain_i.shape (11,6)
       ene_int=ene_i-float(nat0d)/float(nat0b)*E0b-(ene_i_b-E0b)-Efor-edipole # Comment TJ: this is not used; be careful,
       # there is a factor float(nat0d)/float(nat0b) missing in front of (ene_i_b-E0b).

       sigma_tot    =np.vstack([sigma_tot,sigma_i])
       sigma_tot_b  =np.vstack([sigma_tot_b,sigma_i_b])
       strain_tot   =np.vstack([strain_tot,strain_i])
       ene_int_tot  =np.vstack([ene_int_tot,ene_int])
       ene_tot      =np.vstack([ene_tot,ene_i])
       ene_tot_b    =np.vstack([ene_tot_b,ene_i_b])
       ene_tot_corr =np.vstack([ene_tot_corr,ene_corr_i])
  #debug  print 'ene_shape',  ene_tot.shape (77,1), Pij0.shape (1,6), sigma_tot.shpe (77,6)
  if (typerun=='extract') or (typerun=='all'):

    cij=voigt_elastic_matrix_bcc(globalv.valc11,globalv.valc12,globalv.valc44)
    print("C_ij: Elastic matrix in Voigt notation (eV):")
    printMatrixE(cij*V0b/unit_conversion1)

    print("P_ij^0: 0 order of dipole tensor in Voigt notation (eV):")
    Pij0=(vsigma_d0/10.0-vsigma_b0/10.0)*V0b/globalv.unit_conversion1
    printMatrixE(Pij0)

    Efor=formation_energy(nat0b,nat0d,E0b,E0d)
    print("Formation energy  (eV): ")
    print(Efor)

    alpha0 =fit_linalg_leastq (sigma_tot, sigma_tot_b, strain_tot, vsigma_d0, vsigma_b0, V0d, V0b, cij)
    #alpha1 =fit_homemade_leastq (sigma_tot, sigma_tot_b, strain_tot, vsigma_d0, vsigma_b0, V0d, V0b, cij)

    ecorr_sym=func_ecorr_sym (def_type)
    alpha  =fit_scipy_leastq (ene_tot, ene_tot_b, ene_tot_corr, ecorr_sym, strain_tot, Pij0, nat0b, nat0d, E0b, E0d, V0b, V0d, cij, alpha0)

    ii=0
    for def_type in globalv.list_of_deformations:

      test_polar_deformation(based, def_type, def_xp_strain[ii,:], \
                          ene_tot, ene_tot_b, ene_tot_corr, sigma_tot, sigma_tot_b,strain_tot, \
                          E0b, E0d, Efor, sigma0b, sigma0d, nat0b, nat0d,V0b, V0d, \
                          cij,Pij0, \
                          alpha,alpha0)
      ii+=1


  return;



def get_polar_deformation(based,wcname,def_type,a0,E0,sigma0b, sigma0d,xp_vol):
  mode=globalv.mode
  dynamics=globalv.dynamics
  typerun=globalv.typerun
  bulk_deform=globalv.bulk_deform
  ax,ay,az=a0,a0,a0
  limit_sigma=globalv.limit_sigma
  ene_i=np.array([]).reshape(0,1)
  ene_i_b=np.array([]).reshape(0,1)
  ene_corr_i=np.array([]).reshape(0,1)
  sigma_i=np.array([]).reshape(0,6)
  sigma_i_b=np.array([]).reshape(0,6)
  strain_i=np.array([]).reshape(0,6)
  if (mode == 'vasp') :
     setup_param_vasp.init()
  if (mode == 'ndm') :
     setup_param_ndm.init()
  if (mode == 'lammps') :
     setup_param_lammps.init()
  if (mode == 'pwscf') :
     setup_param_pwscf.init()
     prefix=setup_param_pwscf.prefix
#dmcm atoms_per_unit=globalv.atoms_per_unit
#Calculate elastic constants
  #this gives names like deformation001
  cname=wcname+"{0:0>3d}".format(def_type)

  energy = []
  volume = []
  print("... %s"%cname)
  os.system("echo %s > breport.%s"%(cname,mode))
  os.chdir(based)
  if (typerun == 'build'  or ( typerun == 'all') ):

   if globalv.abinitio==0:
     os.system('rm -rf %s'%cname)
     os.mkdir("%s"%cname)

  nvalue=int(10000)
  itest=0
  #debug for i in range(len(xp_vol)):
  #debug   print i, xp_vol[i]
  #debug exit(0)
  for i in range(len(xp_vol)):
    t = xp_vol[i]
    #debug_polar
    #debug_polar print i,t
    #Generate lattice according to strain tensor
    indir=str(nvalue+i)
    indir_b='ebulk_'+str(nvalue+i)
    #One by one version ....
    tempd=based+'/'+cname
    os.chdir(tempd)

    # apply strain
    strain_tensor=apply_strain_woo_nonorm (def_type,t)
    #strain_tensor=apply_strain (def_type,t)
    #the Voight strain was filled ...


    if (typerun == 'build') or (typerun == 'all' ):

      os.mkdir("%s"%str(indir))
      os.chdir("%s"%str(indir))

# building structures files
  # it neeeds: gen_structure, strain_tensor, inpDir, mode, dynamics, ginfile, a0
     #gen_structure([ax,ay,az],strain_tensor,ginFile='structure_dfct.gin')

      write_structure(ax,ay,az, strain_tensor,indir,ginFile='structure_dfct.gin')
      print_strain(xp_vol[0],xp_vol[-1],len(xp_vol), i,xp_vol[i],strain_tensor)

      if bulk_deform=='yes':
         os.chdir('../')

         os.mkdir("%s"%str(indir_b))
         os.chdir("%s"%str(indir_b))

         write_structure(ax,ay,az, strain_tensor,indir_b,ginFile='structure_bulk.gin')
         print_strain(xp_vol[0],xp_vol[-1],len(xp_vol), i,xp_vol[i],strain_tensor)

         os.chdir('../%s'%indir)


      setup_calculation (dirBase=based, dirDeformation=cname, dirLocal=indir)


    if (typerun == 'run' ) or (typerun == 'all') :
      os.chdir(tempd+'/'+indir)
      run_md(mode, tempd, indir,Exceptions=None)

      if bulk_deform=='yes':
        os.chdir(tempd+'/'+indir_b)
        run_md(mode, tempd, indir_b,Exceptions=None)
        os.chdir(tempd+'/'+indir)


    if (typerun == 'extract') or (typerun == 'all') :
      #debug print tempd, indir
      os.chdir(tempd+'/'+indir)
      if ( mode=='ndm') :
        nat=get_natom_from_ndm("out.run")
        ene_local=get_energy_from_ndm("out.run")
        sigma_local=pack_sigma_in_vsigma(get_sigma_from_ndm("out.run"))
        strain_local=np.array(strain_tensor)
        if bulk_deform=='no':
          sigma_local_b=np.array([0,0,0,0,0,0]).reshape(1,6)
          ene_local_b=0.0
        if bulk_deform=='yes':
          os.chdir(tempd+'/'+indir_b)
          ene_local_b = get_energy_from_ndm("out.run")
          sigma_local_b=pack_sigma_in_vsigma(get_sigma_from_ndm("out.run"))
          os.chdir(tempd+'/'+indir)
        if (dynamics=='neb'):
          print("not implemented NEB dynamics with ndm")
          exit(0)


      if ( mode=='lammps') :

        if (dynamics=='neb' and indir[1:5] != 'bulk'):
          nat, ene_local, volume_tmp, sigma_fin = get_all_from_neb_lammps("out.run")
          volume.append(volume_tmp)
        else:
          ene_local = get_energy_from_lammps("out.run",Patched=True)
          volume_tmp = get_volume_from_lammps("out.run",Patched=True)
          volume.append(volume_tmp)
          sigma_fin=get_sigma_from_lammps("out.run",Patched=True)
          nat=get_natom_from_lammps("out.run",Patched=True)

        gen_structure([ax,ay,az],strain_tensor,ginFile='structure_dfct.gin') # in order to recover the beta_matrix
        beta_matrix=globalv.beta_matrix

        sigma_ini=np.matmul(beta_matrix.T,np.matmul(sigma_fin,beta_matrix))
        sigma_local=pack_sigma_in_vsigma(sigma_ini)
        strain_local=np.array(strain_tensor)
        if bulk_deform=='no':
          sigma_local_b=np.array([0,0,0,0,0,0]).reshape(1,6)
          ene_local_b=0.0
        if bulk_deform=='yes':
          os.chdir(tempd+'/'+indir_b)
          ene_local_b = get_energy_from_lammps("out.run",Patched=True)
          sigma_fin_b=get_sigma_from_lammps("out.run",Patched=True)


          sigma_ini_b=np.matmul(beta_matrix.T,np.matmul(sigma_fin_b,beta_matrix))
          sigma_local_b=pack_sigma_in_vsigma(sigma_ini_b)
          os.chdir(tempd+'/'+indir)





      if (mode =='pwscf' ):
        print('polar_deformation for pwscf: this mode not yet implemented')
        exit(0)

      if (mode=='phondy' ):
        print('polar_deformation for phondy: this mode not yet implemented')
        exit(0)

      if ( mode=='vasp' ):
        nat=get_natom_from_vasp_outcar("OUTCAR")
        energy.append(get_energy_from_vasp_outcar("OUTCAR"))
        volume.append(get_volume_from_vasp_outcar("OUTCAR"))
        sigma_local=pack_sigma_in_vsigma(get_sigma_from_vasp("OUTCAR"))

        strain_local=np.array(strain_tensor)
        ene_local=get_energy_from_vasp_outcar("OUTCAR")


        if bulk_deform=='no':
          sigma_local_b=np.array([0,0,0,0,0,0]).reshape(1,6)
          ene_local_b=0.0
        if bulk_deform=='yes':
          os.chdir(tempd+'/'+indir_b)
          #debug os.system('pwd')
          ene_local_b=get_energy_from_vasp_outcar("OUTCAR")
          sigma_local_b=pack_sigma_in_vsigma(get_sigma_from_vasp("OUTCAR"))
          os.chdir(tempd+'/'+indir)

      ecorr=0.0
      if (globalv.aneto=='yes') :
           if bulk_deform=='no':
             ecorr = pack_aneto (dirbased=based, bulkdir='0bulk', sigma_bulk=sigma0b, defectdir=cname+'/'+indir,sigma_defect=sigma_local)
           if bulk_deform=='yes':
             ecorr = pack_aneto (dirbased=based, bulkdir=cname+'/'+indir_b , sigma_bulk=sigma_local_b, defectdir=cname+'/'+indir,sigma_defect=sigma_local)
      else:
           ecorr=0.0




    #we pack all the sigmas ...
      #if (i==0):
      #     sigma_i=sigma_local
      #     strain_i=strain_local
      #if (i>0):
      sigma_i=np.vstack([sigma_i,sigma_local])
      sigma_i_b=np.vstack([sigma_i_b,sigma_local_b])
      strain_i=np.vstack([strain_i,strain_local])
      ene_i=np.vstack([ene_i,ene_local])
      ene_i_b=np.vstack([ene_i_b,ene_local_b])
      ene_corr_i=np.vstack([ene_corr_i,ecorr])

  return ene_i, ene_i_b, ene_corr_i, sigma_i,sigma_i_b,strain_i




def  test_polar_deformation(based, def_type, xp_strain, \
                          ene_i, ene_i_b, ene_i_corr, sigma_i, sigma_i_b,strain_i, \
                          E0b, E0d, Efor, sigma0b, sigma0d, nat0b, nat0d,V0b, V0d, \
                          cij,Pij0, \
                          alpha,alpha_sigma):

  plotlincorr = True # TJ: to remove the dipolar contribution
  plotpercent = True
  if plotpercent:
    factpercent = 100.
  else:
    factpercent = 1.
  
  mode=globalv.mode
  dynamics=globalv.dynamics
  typerun=globalv.typerun
  bulk_deform=globalv.bulk_deform

  os.chdir(based)
  nvalue=int(10000)
  itest=0

  cname="deformation"+"{0:0>3d}".format(def_type)

  idef=0
  for i_def_type in globalv.list_of_deformations:
    if (i_def_type==def_type):
      local_ene_i=ene_i[idef:idef+globalv.nt_strain,:]
      local_ene_i_corr=ene_i_corr[idef:idef+globalv.nt_strain,:]
      local_sigma_i=sigma_i[idef:idef+globalv.nt_strain,:]
      local_strain_i=strain_i[idef:idef+globalv.nt_strain,:]
      local_sigma_i_b=sigma_i_b[idef:idef+globalv.nt_strain,:]
      local_ene_i_b=ene_i_b[idef:idef+globalv.nt_strain,:]
    idef+=globalv.nt_strain


  unit_conversion1 = 160.217646


  al=np.matmul(local_strain_i,np.matmul(cij.T,local_strain_i.T))*V0b/(2.0*unit_conversion1)

  # the elastic bulk energy is accounted only from elastic constants ...
  y_atomic=local_ene_i[:,0]-float(nat0d)/float(nat0b)*E0b-Efor-np.diagonal(al)

  imiddle=int((globalv.nt_strain-1)/2)

  ecorr_sym=func_ecorr_sym (def_type)


  if bulk_deform=='yes':
     #the elastic energy of bulk deformation is computed atomisticaly ene_i_b - E0b
     y_atomic_full=local_ene_i[:,0]-E0d \
                   -float(nat0d)/float(nat0b)*(local_ene_i_b[:,0]-E0b)+ecorr_sym
     #the elastic energy of bulk deformation is computed atomisticaly ene_i_b - E0b + aneto correction
     y_atomic_full_corr=local_ene_i[:,0]-E0d-local_ene_i_corr[:,0]+local_ene_i_corr[imiddle,0] \
                   -float(nat0d)/float(nat0b)*(local_ene_i_b[:,0]-E0b)+ecorr_sym
              #     -1.0*(local_ene_i_b[:,0]-E0b)
     y_atomic_full_corr[imiddle]=y_atomic_full_corr[imiddle]-ecorr_sym
     y_atomic_full[imiddle]=y_atomic_full[imiddle]-ecorr_sym
     #debug
     y_atomic_only_bulk=-float(nat0d)/float(nat0b)*(local_ene_i_b[:,0]-E0b)
     y_atomic_only_dfct=-(local_ene_i[:,0]-E0d)



  al0=np.matmul(-Pij0,local_strain_i.T)
  y_atomic_full_corr=y_atomic_full_corr
  if plotlincorr:
    y_atomic_full_corr = y_atomic_full_corr-al0[0,:]


  #al0=al
  #al0.shape 1,1 , Pij0.shape 1,6, strain_i.shape 11,6
  aint=np.matmul(local_strain_i,np.matmul(-alpha,local_strain_i.T))/2.0
  aint_sigma=np.matmul(local_strain_i,np.matmul(-alpha_sigma,local_strain_i.T))/2.0
  al1=np.matmul(-Pij0,local_strain_i.T)+np.diagonal(aint)
  al1_sigma=np.matmul(-Pij0,local_strain_i.T)+np.diagonal(aint_sigma)

  #print Pij0
  y_elastic0=al0[0,:]
  y_elastic1=al1[0,:]
  if plotlincorr:
    y_elastic1 = y_elastic1 - al0[0,:]
  y_elastic1_sigma=al1_sigma[0,:]
  if plotlincorr:
    y_elastic1_sigma = y_elastic1_sigma - al0[0,:]

  #debug print y_elastic0.shape, y_elastic1.shape


  emin=xp_strain[0]-0.01/100.0
  emax=xp_strain[-1]+0.01/100.0
  xepsilon=np.linspace(emin,emax,num=20,endpoint=True)
  #y_atomic_smooth=interp1d(xp_vol/100.0,y_atomic,kind='cubic')
#elastic  y_atomic_smooth=CubicSpline(xp_strain,y_atomic)
  if bulk_deform=='yes':
    y_atomic_full_smooth=CubicSpline(xp_strain,y_atomic_full)
    y_atomic_full_corr_smooth=CubicSpline(xp_strain,y_atomic_full_corr)

    y_atomic_only_bulk_smooth=CubicSpline(xp_strain,y_atomic_only_bulk)
    y_atomic_only_dfct_smooth=CubicSpline(xp_strain,y_atomic_only_dfct)

  y_elastic0_smooth=CubicSpline(xp_strain,y_elastic0)
  y_elastic1_smooth=CubicSpline(xp_strain,y_elastic1)
  y_elastic1_sigma_smooth=CubicSpline(xp_strain,y_elastic1_sigma)
  fig=pl.figure()
  ax=fig.add_subplot(111)
#elastic  ax.plot(xp_strain,y_atomic,'o',markersize=6,color='slateblue',label='atomic-elastic')

  pl.rc('legend',**{'fontsize':20})
  if bulk_deform=='yes':
    #debug ax.plot(xp_strain,y_atomic_full,'o',markersize=6,color='green',label='atomic-full')
    #ax.plot(xepsilon,y_atomic_full_smooth(xepsilon),linestyle='-',markersize=6,color='green')

    ##debug
    #ax.plot(xp_strain,y_atomic_only_bulk,'o',markersize=6,color='skyblue',label='atomic-bulk')
    #ax.plot(xepsilon,y_atomic_only_bulk_smooth(xepsilon),'-',markersize=6,color='skyblue')
    ##debug
    #ax.plot(xp_strain,y_atomic_only_dfct,'o',markersize=6,color='m',label='atomic-dfct')
    #ax.plot(xepsilon,y_atomic_only_dfct_smooth(xepsilon),'-',markersize=6,color='m')


    ax.plot(factpercent*xp_strain,y_atomic_full_corr,'o',markersize=6,color='lightgreen',label='DFT')
    ax.plot(factpercent*xepsilon,y_atomic_full_corr_smooth(xepsilon),linestyle='-',markersize=6,color='lightgreen')
#elastic   ax.plot(xepsilon,y_atomic_smooth(xepsilon),linestyle='-',markersize=6,color='slateblue')

  if not plotlincorr:
    ax.plot(factpercent*xp_strain,y_elastic0,'v',markersize=6,color='tomato',label=r"$P_{ij} \approx P_{ij}^{(0)}$")
    ax.plot(factpercent*xepsilon,y_elastic0_smooth(xepsilon),linestyle='dotted',markersize=6,color='tomato')

  ax.plot(factpercent*xp_strain,y_elastic1,'s',markersize=6,color='plum', label=r"$P_{ij} \approx P_{ij}^{(0)} +  \alpha_{ijkl}^E\epsilon_{kl}$")
  ax.plot(factpercent*xepsilon,y_elastic1_smooth(xepsilon),linestyle='--',markersize=6,color='plum')

  ax.plot(factpercent*xp_strain,y_elastic1_sigma,'^',markersize=6,color='skyblue', label=r"$P_{ij} \approx P_{ij}^{(0)} +  \alpha_{ijkl}^{\sigma}\epsilon_{kl}$")
  ax.plot(factpercent*factpercent*xepsilon,y_elastic1_sigma_smooth(xepsilon),linestyle='--',markersize=6,color='skyblue')



  #ax.set_ylabel(r" $E^{b+d}(\epsilon) - E^{b+d}(\epsilon=0) - E^{b}(\epsilon) + E^b(\epsilon=0) $ (eV)", color='slateblue', fontsize='20')
  ax.set_ylabel(r" $ \Delta E (\epsilon) $ (eV)", fontsize='20')
  if plotpercent:
    ax.set_xlabel(r"$\epsilon$ (\%)", fontsize='20')
  else:
    ax.set_xlabel(r"$\epsilon$", fontsize='20')
  ax.set_xlim([factpercent*emin,factpercent*emax])
  ax.tick_params(axis='x',labelsize='16')
  ax.tick_params(axis='y',labelsize='16')
  ax.legend(loc='best')
  os.chdir(based)
  defsize=fig.get_size_inches()
  fig.set_size_inches( (defsize[0]*1.5,defsize[1]*1.5) )
  fig.savefig('%s.png'%cname)
  #fig.savefig('%s.eps'%cname)
  #plot_show
#  pl.show()

  return





def  test_vasp_polar(based,wcname):

 mode=globalv.mode
 typerun=globalv.typerun
 setup_param_vasp.init()
 prefix=setup_param_vasp.prefix
 nt_strain=globalv.nt_strain
 bulk_deform=globalv.bulk_deform
 #testing pre files

 nvalue=int(10000)
 for def_type in globalv.list_of_deformations:
   cname=wcname+"{0:0>3d}".format(def_type)

   for i in range(nt_strain):
     indir=str(nvalue+i)
     #One by one version ....
     tempd=based+'/'+cname
     #debug print tempd, based, cname
     os.chdir(tempd)
     os.chdir("%s"%str(indir))
     if os.path.exists('OUTCAR'):
       os.system('pwd')
       iout = test_if_vasp_terminated_correctly("OUTCAR")
       if iout == 0:
         print('OUTCAR unfinished: %s'%(tempd+'/'+str(indir)))
     else:
         print("OUTCAR missing     : %s"%(tempd+'/'+str(indir)))


   if (bulk_deform=='yes'):
     for i in range(nt_strain):
       indir=str(nvalue+i)
       #One by one version ....
       tempd=based+'/'+cname
       #debug print tempd, based, cname
       os.chdir(tempd)
       os.chdir("%s"%('ebulk_'+str(indir)))

       if os.path.exists('OUTCAR'):
         iout = test_if_vasp_terminated_correctly("OUTCAR")
         if iout == 0:
           print('OUTCAR unfinished: %s'%(tempd+'/'+'ebulk_'+str(indir)))
       else:
           print("OUTCAR missing     : %s"%(tempd+'/'+'ebulk_'+str(indir)))

 exit(0)
 return ;

def  test_lammps_polar(based,wcname):

 mode=globalv.mode
 typerun=globalv.typerun
 setup_param_lammps.init()
 nt_strain=globalv.nt_strain
 #testing pre files

 nvalue=int(10000)
 for def_type in globalv.list_of_deformations:
   cname=wcname+"{0:0>3d}".format(def_type)

   for i in range(nt_strain):
     indir=str(nvalue+i)
     #One by one version ....
     tempd=based+'/'+cname
     #debug print tempd, based, cname
     os.chdir(tempd)
     os.chdir("%s"%str(indir))

     if os.path.exists('out.run'):
       iout = test_if_lammps_terminated_correctly("out.run")
       if iout == 0:
         print('our.run unfinished  : %s'%(tempd+'/'+str(indir)))
         run_lammps(tempd, indir,"out.run")
     else:
         print('out.run missing     : %s'%(tempd+'/'+str(indir)))
         run_lammps(tempd, indir,"out.run")
 exit(0)
 return ;


def apply_strain (def_type,t):

    if def_type==1 :
         strain_tensor=[t,t,t,0,0,0]
    elif def_type==2 :
         strain_tensor=[-2*t,t,t,0,0,0]
    elif def_type==3 :
         strain_tensor=[t,0,-t,2.0*t,0,0]
    elif def_type==4 :
         strain_tensor=[0,0,0,2.0*t,0,2.0*t]
    elif def_type==5 :
         strain_tensor=[0,0,0,0,2.0*t,0]
    elif def_type==6 :
         strain_tensor=[t,-t,1.0/(1.0-t**2)-1,0,0,0]
    elif def_type==7 :
         strain_tensor=[0,0,t**2/(4.0-t**2),0,0,2.0*t]
    elif def_type==8 :
         strain_tensor=[0,0,0,0,0,2.0*t]
    elif def_type==9 :
         strain_tensor=[0,t,0,0,0,0]
    elif def_type==10 :
         strain_tensor=[0,0,t,0,0,0]
    elif def_type==11 :
         strain_tensor=[0,0,0,0,2.0*t,2.0*t]
    elif def_type==12 :
         strain_tensor=[0,t,0,2.0*t,2.0*t,-2.0*t]
    else:
           print("Deformation type not known. def_type=%s"%def_type)

    return strain_tensor;

def apply_strain_no_two (def_type,t):

    if def_type==1 :
         strain_tensor=[t,t,t,0,0,0]
    elif def_type==2 :
         strain_tensor=[-2*t,t,t,0,0,0]
    elif def_type==3 :
         strain_tensor=[t,0,-t,t,0,0]
    elif def_type==4 :
         strain_tensor=[0,0,0,t,0,t]
    elif def_type==5 :
         strain_tensor=[0,0,0,0,t,0]
    elif def_type==6 :
         strain_tensor=[t,-t,1.0/(1.0-t**2)-1,0,0,0]
    elif def_type==7 :
         strain_tensor=[0,0,t**2/(4.0-t**2),0,0,t]
    elif def_type==8 :
         strain_tensor=[0,0,0,0,0,t]
    elif def_type==9 :
         strain_tensor=[0,t,0,0,0,0]
    elif def_type==10 :
         strain_tensor=[0,0,t,0,0,0]
    elif def_type==11 :
         strain_tensor=[0,0,0,0,t,t]
    elif def_type==12 :
         strain_tensor=[0,t,0,t,t,-t]
    else:
           print("Deformation type not known. def_type=%s"%def_type)

    return strain_tensor;


def apply_strain_woo (def_type,t):

    if def_type==1 :
         strain_tensor=[t,t,t,0,0,0]/np.sqrt(3.0)
    elif def_type==2 :
         strain_tensor=[t,-t,0,0,0,0]/np.sqrt(2.0)
    elif def_type==3 :
         strain_tensor=[t,t,-2.0*t,0,0,0]/np.sqrt(6.0)
    elif def_type==4 :
         strain_tensor=[0,0,0,2.0*t,0,0]/np.sqrt(2.0)
    elif def_type==5 :
         strain_tensor=[0,0,0,0,2.0*t,0]/np.sqrt(2.0)
    elif def_type==6 :
         strain_tensor=[0.,0.,0.,0,0,2.0*t]/np.sqrt(2.0)
    elif def_type==7 :
         strain_tensor=[0,0,t**2/(4.0-t**2),0,0,2.0*t]
    elif def_type==8 :
         strain_tensor=[t,-t,1.0/(1.0-t**2)-1,0,0,0]
    elif def_type==9 :
         strain_tensor=[0.,0.,0.,0.,2.0*t,2.0*t]/np.sqrt(4.0)
    elif def_type==10 :
         strain_tensor=[0.,0.,0.,0.,2.0*t,-2.0*t]/np.sqrt(4.0)
    elif def_type==11 :
         strain_tensor=[0,t,0,2.0*t,2.0*t,-2.0*t]/np.sqrt(5.0)
    else:
           print("Deformation type not known. def_type=%s"%def_type)

    return strain_tensor;


def apply_strain_woo_nonorm (def_type,t):

    if def_type==1 :
         strain_tensor=[t,t,t,0,0,0]#/np.sqrt(3.0)
    elif def_type==2 :
         strain_tensor=[t,-t,0,0,0,0]#/np.sqrt(2.0)
    elif def_type==3 :
         strain_tensor=[t,t,-2.0*t,0,0,0]#/np.sqrt(6.0)
    elif def_type==4 :
         strain_tensor=[0,0,0,2.0*t,0,0]#/np.sqrt(2.0)
    #old elif def_type==5 :
    #old     strain_tensor=[0,0,0,0,2.0*t,0]#/np.sqrt(2.0)
    elif def_type==5 :
         strain_tensor=[0,0,0,t,t,0]#/np.sqrt(2.0)
    elif def_type==6 :
         strain_tensor=[0.,0.,0.,0,0,2.0*t]#/np.sqrt(2.0)
    #old elif def_type==7 :
    #old      strain_tensor=[0,0,t**2/(4.0-t**2),0,0,2.0*t]
    elif def_type==7 :
         #strain_tensor=[t,0,0,0,0,t]
         strain_tensor=[t,0,0,0,0,2*t]
         # What Thomas suggests strain_tensor=[t,0,0,0,0,2.0*t]
    #old elif def_type==8 :
    #old     strain_tensor=[t,-t,1.0/(1.0-t**2)-1,0,0,0]
    elif def_type==8 :
          #strain_tensor=[0,0,t,0,0,t]
         strain_tensor=[0,0,t,0,0,2*t]
         # what Thomas suggests strain_tensor=[0,0,t,0,0,2.0*t]
    #elif def_type==9 :
    #     strain_tensor=[0.,0.,0.,0.,2.0*t,2.0*t]#/np.sqrt(4.0)
    elif def_type==9 :
         strain_tensor=[t,t,-t,0.,0.,0.]#/np.sqrt(4.0)
    elif def_type==10 :
         strain_tensor=[0.,0.,0.,0.,2.0*t,-2.0*t]#/np.sqrt(4.0)
    elif def_type==11 :
         strain_tensor=[0,t,2.0*t,0,0,0]#/np.sqrt(5.0)
    elif def_type==23 :
         strain_tensor=[0.0,t,-t,0,0,0]#/np.sqrt(6.0)



    else:
           print("Deformation type not known. def_type=%s"%def_type)

    return strain_tensor;




# Function : i_use_line
# author : Denise Carpentier
# date : 23/11/17
# description : i_use_line identifies if the line should be rejected,
# used for training or used for validation.
# input : i_line is an integer between 0 and 76 (function built for 7 strains sampled with 11 points)
# output : i_use is and integr : i_use = i_reject => line rejected.
#                                i_use = i_train => line used for training
#                                i_use = i_valid => line used for validation
def i_use_line(i_line,nlen):
  # define the wanted values :
  i_reject = 0 # => line rejected.
  i_train = 1  # => line used for training
  i_valid = 2  # => line used for validation
  number_of_strains = 7
  number_of_samples = nlen

  if i_line < 0 or i_line >= (number_of_strains*number_of_samples) :
    print('error! the function i_use_line was intended to treat line indeces between 0 and'+str((number_of_strains*number_of_samples)-1)+' (function built for '+str(number_of_strains)+' strains sampled with '+str(number_of_samples)+' points)! Index given :'+str(i_line))
    exit
  else:
    i_strain = i_line/number_of_samples
    # print 'Line '+str(i_line)+' corresponds to strain type '+str(i_strain)
    i_sample = i_line%number_of_samples
    # print 'it is the sample '+str(i_sample)
    # example for 11 samples :
    #             | sample 0
    #             | sample 1  <-- validation
    #             | sample 2
    #             | sample 3
    #             | sample 4  <-- rejected
    # strain type | sample 5
    #             | sample 6  <-- rejected
    #             | sample 7
    #             | sample 8
    #             | sample 9  <-- validation
    #             | sample 10
    if i_sample == 1 or i_sample == number_of_samples-2 :
      i_use = i_valid
    elif i_sample == number_of_samples/2 - 1 or i_sample == number_of_samples/2 + 1:
      i_use = i_reject
    else:
      i_use = i_train

  return i_use
