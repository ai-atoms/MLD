
def  test_polar_deformation(based,wcname,def_type,a0,E0,xp_vol,step_xp_vol ,\
                          cij,\
                          V0b,\
                          E0b,\
                          Efor, \
                          nat0b,\
                          nat0d,\
                          Pij0, \
                          alpha):

  mode=globalv.mode
  dynamics=globalv.dynamics
  typerun=globalv.typerun
  bulk_deform=globalv.bulk_deform
  ax,ay,az=a0,a0,a0 
  limit_sigma=globalv.limit_sigma
  sigma_i=np.array([]).reshape(0,6)
  sigma_i_p=np.array([]).reshape(0,6)
  strain_i=np.array([]).reshape(0,6)
  cname=wcname+"{0:0>3d}".format(def_type)

  os.chdir(based)
  energy=[]
  volume=[]
  nvalue=int(10000)
  itest=0
  y_atomic=[]
  y_atomic_full=[]
  y_elastic0=[]
  y_elastic1=[]
  for i in range(len(xp_vol)):
    t = xp_vol[i]
    #debug_polar print 'e,',i,t
    #Generate lattice according to strain tensor
    indir=str(nvalue+i)
    indir_b='ebulk_'+str(nvalue+i)
    #One by one version .... 
    tempd=based+'/'+cname
    #print "tmpd  %s"%tempd
    os.chdir(tempd)
     
    # apply strain
    strain_tensor=apply_strain_woo_nonorm (def_type,t)
    #strain_tensor=apply_strain (def_type,t)
    #the Voight strain was filled ...
    
    if (typerun == 'extract') or (typerun == 'all') :
      #debug print tempd, indir 
      os.chdir(tempd+'/'+indir)
      if (mode=='ndm') :
       nat=get_natom_from_ndm("out.run")
       ene_i=get_energy_from_ndm("out.run")
       sigma_local=pack_sigma_in_vsigma(get_sigma_from_ndm("out.run"))
       strain_local=np.array(strain_tensor)
      if ( mode=='lammps') :

        if (dynamics=='neb' and indir[1:5] != 'bulk'):
          nat, ene_i, volume_tmp, sigma_fin = get_all_from_neb_lammps("out.run")
          volume.append(volume_tmp)
        else:
          ene_i = get_energy_from_lammps("out.run",Patched=True)
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

        if bulk_deform=='yes':
          os.chdir(tempd+'/'+indir_b)
          ene_i_b = get_energy_from_lammps("out.run",Patched=True)
          sigma_fin_b=get_sigma_from_lammps("out.run",Patched=True)
 
 
          sigma_ini_b=np.matmul(beta_matrix.T,np.matmul(sigma_fin_b,beta_matrix))
          sigma_local_b=pack_sigma_in_vsigma(sigma_ini_b)
          os.chdir(tempd+'/'+indir)



      if ( mode=='pwscf') :
       nat=get_natom_from_pwscf_outcar("out.run") 
       energy.append(get_energy_from_pwscf_outcar("out.run")/float(nat))
       volume.append(get_volume_from_pwscf_outcar("out.run")/float(nat))
       get_sigma_from_pwscf("out.run")
      if (mode == 'bphondy' or mode=='phondy' ):
       print 'polar mode, phondy: this mode not yet implemented'
       exit(0)
      if ( mode=='vasp' ):
  
       nat=get_natom_from_vasp_outcar("OUTCAR")
       ene_i=get_energy_from_vasp_outcar("OUTCAR")
       #volume.append(get_volume_from_vasp_outcar("OUTCAR"))
       sigma_local=pack_sigma_in_vsigma(get_sigma_from_vasp("OUTCAR"))
       strain_local=np.array(strain_tensor)

       if bulk_deform=='no':
         sigma_local_b=np.array([0,0,0,0,0,0]).reshape(1,6)

       if bulk_deform=='yes':
         os.chdir(tempd+'/'+indir_b)
         ene_i_b = get_energy_from_vasp_outcar("OUTCAR")
         sigma_ini_b=get_sigma_from_vasp("OUTCAR")
         sigma_local_b=pack_sigma_in_vsigma(sigma_ini_b)
         os.chdir(tempd+'/'+indir)



      unit_conversion1 = 160.217646
      strain_i=np.array([strain_local])
      #print np.shape(Pij0) 
      al=np.matmul(strain_i,np.matmul(cij.T,strain_i.T))*V0b/(2.0*unit_conversion1)

      y_atomic.append(ene_i-float(nat0d)/float(nat0b)*E0b-Efor-al[0,0])
      if bulk_deform=='yes':
        y_atomic_full.append(ene_i-float(nat0d)/float(nat0b)*E0b-Efor-float(nat0d)/float(nat0b)*(ene_i_b-E0b))

        #y_atomic_full.append((ene_i_b-E0b)-al[0,0])
        #y_atomic_full.append((ene_i_b-E0b))
  
   
      al0=np.matmul(-Pij0,strain_i.T)
      #al0=al
      #al0.shape 1,1 , Pij0.shape 1,6, strain_i.shape 1,6
      al1=np.matmul(-Pij0,strain_i.T)-np.matmul(strain_i,np.matmul(alpha,strain_i.T))/2.0

      #print Pij0


 
      y_elastic0.append(al0[0,0])
      y_elastic1.append(al1[0,0])
      
      #debug print Efor+np.matmul(strain_i,np.matmul(cij.T,strain_i.T))*V0b/(2.0*unit_conversion1)
      

  emin=xp_vol[0]-0.01/100.0
  emax=xp_vol[-1]+0.01/100.0
  xepsilon=np.linspace(emin,emax,num=20,endpoint=True)
  #y_atomic_smooth=interp1d(xp_vol/100.0,y_atomic,kind='cubic')
  y_atomic_smooth=CubicSpline(xp_vol,y_atomic)
  if bulk_deform=='yes':
    y_atomic_full_smooth=CubicSpline(xp_vol,y_atomic_full)
  y_elastic0_smooth=CubicSpline(xp_vol,y_elastic0)
  y_elastic1_smooth=CubicSpline(xp_vol,y_elastic1)
  fig=pl.figure()
  ax=fig.add_subplot(111)
  ax.plot(xp_vol,y_atomic,'o',markersize=6,color='slateblue',label='atomic-elastic') 
  if bulk_deform=='yes':
    ax.plot(xp_vol,y_atomic_full,'o',markersize=6,color='green',label='atomic-full') 
    ax.plot(xepsilon,y_atomic_full_smooth(xepsilon),linestyle='-',markersize=6,color='green') 
  ax.plot(xepsilon,y_atomic_smooth(xepsilon),linestyle='-',markersize=6,color='slateblue') 
  ax.plot(xp_vol,y_elastic0,'v',markersize=6,color='tomato',label=r"$P_{ij}^{(0)}$") 
  ax.plot(xepsilon,y_elastic0_smooth(xepsilon),linestyle='dotted',markersize=6,color='tomato') 
  ax.plot(xp_vol,y_elastic1,'s',markersize=6,color='plum', label=r"$P_{ij}^{(1)}$") 
  ax.plot(xepsilon,y_elastic1_smooth(xepsilon),linestyle='--',markersize=6,color='plum') 
  ax.set_ylabel(r" $E^{b+d}(\epsilon) - E^b(0)-E_{for}-\frac{1}{2}C_{ijkl}\epsilon_{ij}\epsilon_{kl}$ (eV)", color='slateblue',fontsize='20')
  ax.set_xlabel(r"$\epsilon$ (strain)")
  ax.set_xlim([emin,emax])
  ax.tick_params(axis='x',labelsize='16')
  ax.tick_params(axis='y',labelsize='16')
  ax.legend(loc='best') 
  os.chdir(based)
  defsize=fig.get_size_inches()
  fig.set_size_inches( (defsize[0]*1.5,defsize[1]*1.5) )
  fig.savefig('%s.png'%cname)
  fig.savefig('%s.eps'%cname)
  #plot_show 
  pl.show()

  return 



#Denise 


# Function : v_pack_e
# author : Denise Carpentier
# date : 23/11/17
# description :v_pack_e transforms an array [e0, e1, ..., e5] in an array
# [v0, v1, ..., v20]  where v0 = e0^2, v1 = e0e1, ..., v6 = e1^2, ..., v10 = e1e5, ..., v20 = e5^2.
# input : strain_line is a 1-D array of length 6
# output : v_strain_tot is a 1-D array of length 21.
def v_pack_e(strain_line):
  if len(np.shape(strain_line)) > 1 or len(strain_line) != 6:
    print 'error! the function v_pack_e was intended to treat a strain tensor in voigt notation (an array of length 6) ! Given array of shape '+str(np.shape(strain_line))
  else:
    v_strain_tot = np.zeros(21)
    kk = 0
    for ii in range(6):
      for jj in range(ii,6):
        v_strain_tot[kk] = strain_line[ii]*strain_line[jj]
        #print 'v_['+str(kk)+'] = e_'+str(ii)+'*e_'+str(jj)
        kk = kk+1

  return v_strain_tot 



# author : Denise Carpentier
# date : 23/11/17
# description : a_pack_beta converts a 1-D array of length 21
# [a00, a01, ..., a05, a11, ..., a15, ..., a55] to a array of shape (6,6)
# using symetry [[a00, a01, .., a05],[a01, a11, a12, ..., a15], ..., [a05, a15, ..., a55]].
# input : beta_line is a 1-D array of length 21
# output : alpha is a 2D-array of shape (6,6)
def a_pack_beta(beta_line):

  if len(beta_line) != 21 :
    print 'error! the function a_pack_beta is intended to unpack an array of length 21! Given array of shape '+str(np.shape(beta_line))
    exit
  else :
      
    alpha = np.zeros([6,6])

    kk=0
      
    for ii in range(6):
      for jj in range(6):
        if jj < ii :
          alpha[ii,jj] = alpha[ii,jj]
        else:
          alpha[ii,jj] = beta_line[kk]
          kk = kk+1


    print 'Unpacked polarizability matrix :'
    printMatrixE(alpha)
  return alpha

 

def old_func():



    if globalv.bulk_deform=='yes':

      print 'From here, I work:'
      print 'ene_tot.shape, V0b,strain_tot[0].shape'
      print ene_tot.shape, V0b,strain_tot.shape

      #alpha = (-xsol) *V0b/unit_conversion1
            
      print beta0
            
      def function_objective(beta,x):
        E_obj = 0.
        for ia in range(len(beta)):
          E_obj = E_obj+0.5*V0b*beta[ia]*x[ia]
        return E_obj
          
      def func_beta_error(beta, y, x):
        Error = y - function_objective(beta,x[0])
        return Error

      energy_used = []
      strain_line = []
      for ii in range(globalv.nt_strain,len(ene_int_tot)):
        if (i_use_line(ii,globalv.nt_strain) != 0):
          energy_used.append(ene_int_tot[ii])
          v_line = v_pack_e(strain_tot[ii])
          strain_line.append(v_line)
        else:
          print('Rejected '+str(ii))

      E_diff = func_beta_error(beta0,ene_int_tot,np.array(strain_line))
      
      print 'Initial guess energy', E_diff, np.shape(E_diff)

      print np.shape(ene_int_tot[0]), np.shape(strain_line), np.shape(beta0)
      
      beta_new = sp.optimize.least_squares(func_beta_error, np.array(beta0), method='trf', args = (np.array(energy_used[0]),np.array(strain_line)))
      print 'From what I tried, I got:'
      
      print(beta_new.x)
      alpha_energy=a_pack_beta(beta_new.x*V0b/unit_conversion1)

      E_diff = func_beta_error(beta_new.x,ene_int_tot,np.array(strain_line))
      
      print 'Final value energy', E_diff
# end Denise ---------------------------------------------------------                                               


def fit_homemade_validation_leastq (sigma_tot, sigma_tot_b, strain_tot, vsigma_d0, vsigma_b0, V0d, V0b, cij):


#    
    at=np.matmul(cij.T,strain_tot.T)

    #Please pay attention to the sign of sigma in the output file.

    ymat_e=-(sigma_tot-2.0*vsigma_d0+2.0*vsigma_b0)/10.0-at.T
    if globalv.bulk_deform=='yes':
      #Please pay attention to the sign of sigma in the output file.
      # for ndm sigma_elastic=-sigma_output 
      ymat=-((sigma_tot-sigma_tot_b)-2.0*vsigma_d0+2.0*vsigma_b0)/10.0


    if globalv.bulk_deform=='yes':

      x_e_train=np.array([]).reshape(0,6)
      x_e_valid=np.array([]).reshape(0,6)
      y_train=np.array([]).reshape(0,6)
      y_valid=np.array([]).reshape(0,6)
      for i in range(len(strain_tot)):
            if i_use_line(i,nlen)==1:
               x_e_train=np.vstack([x_e_train, strain_tot[i]])
               y_train=np.vstack([y_train,ymat[i]])
            if i_use_line(i,nlen)==2:           
               x_e_valid=np.vstack([x_e_valid, strain_tot[i]])
               y_valid=np.vstack([y_valid,ymat[i]])


      print x_e_train.shape, y_train.shape
      lambda_p=0.0000
      wp_lambda=  np.matmul(np.matmul(np.linalg.inv(np.matmul(x_e_train.T,x_e_train)+lambda_p*np.identity(6)),x_e_train.T),y_train)
      print 'OTHER NEW VERSION'
      betamm = (-wp_lambda) *V0b/unit_conversion1
      printMatrixE(betamm)



#      #regularization fit in 21 space:
      
#      x_e_train=np.array([]).reshape(0,21)
#      x_e_valid=np.array([]).reshape(0,21)
#      y_train=np.array([]).reshape(0,1)
#      y_valid=np.array([]).reshape(0,1)
#      for i in range(len(strain_tot)):
#            if i_use_line(i,nlen)==1:
#               pack_strain_train=v_pack_e(strain_tot[i])
#               x_e_train=np.vstack([x_e_train, pack_strain_train])
#               y_train=np.vstack([y_train,ene_tot[i]])
#            if i_use_line(i,nlen)==2:           
#               pack_strain_validate=v_pack_e(strain_tot[i])
#               x_e_valid=np.vstack([x_e_valid, pack_strain_validate])
#               y_valid=np.vstack([y_valid,ene_tot[i]])
#
#     
#      #print y_train.shape (49,1) , x_e_train.shape (49,21) 
#      wp=  np.matmul(np.matmul(np.linalg.inv(np.matmul(x_e_train.T,x_e_train)),x_e_train.T),y_train)
#      print w, wp.shape

      #pack for training:
     
     
      #print "alpha_ijkl without C_ijkl and with weights : polarizibility matrix in Voigt notation (eV):"
      #alpha = (-xsolw) *V0b/unit_conversion1
      #printMatrixE((-xsolw) *V0b/unit_conversion1)
      #testing the energy ... 
 

