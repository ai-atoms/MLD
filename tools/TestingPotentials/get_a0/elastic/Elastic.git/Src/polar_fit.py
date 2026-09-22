import numpy as np
import globalv
from tools import *
from scipy.optimize import least_squares


def fit_linalg_leastq (sigma_tot, sigma_tot_b, strain_tot, vsigma_d0, vsigma_b0, V0d, V0b, cij):

#
    at=np.matmul(cij.T,strain_tot.T)

    #Please pay attention to the sign of sigma in the output file.
    # for ndm sigma_elastic=-sigma_output
    ymat=-(sigma_tot-2.0*vsigma_d0+2.0*vsigma_b0)/10.0-at.T

    #debug print 'ymat', np.shape(ymat),np.shape(strain_tot.T)
    xsol, residue, ranka, sa = np.linalg.lstsq( strain_tot,ymat, rcond=None)


    #print PrettyTable(xsol)
    #print DataFrame(xsol*V0b/unit_conversion1)
    #Please pay attention to the sign of sigma in the output file.
    # for ndm sigma_elastic=-sigma_output

    print("alpha_ijkl using C_ijkl: polarizibility matrix in Voigt notation (eV):")
    alpha0=-xsol*V0b/globalv.unit_conversion1
    printMatrixE(-xsol*V0b/globalv.unit_conversion1 )

    if globalv.bulk_deform=='yes':
      #Please pay attention to the sign of sigma in the output file.
      # for ndm sigma_elastic=-sigma_output
      ymat=-((sigma_tot-sigma_tot_b)-2.0*vsigma_d0+2.0*vsigma_b0)/10.0

      # Take weigths into account in the lstsq function:

      xsol, residue, ranka, sa = np.linalg.lstsq(strain_tot , ymat, rcond=None)

      #debug printMatrixE(xsol)
      #printMatrixE((cij-xsol) *V0b/unit_conversion1)
      print("alpha_ijkl without C_ijkl: polarizibility matrix in Voigt notation (eV):")
      alpha = (-xsol) *V0b/globalv.unit_conversion1
      printMatrixE((-xsol) *V0b/globalv.unit_conversion1)


    if globalv.bulk_deform=='yes':
      return  alpha
    else:
      return  alpha0


def fit_homemade_leastq (sigma_tot, sigma_tot_b, strain_tot, vsigma_d0, vsigma_b0, V0d, V0b, cij):

#
    at=np.matmul(cij.T,strain_tot.T)

    #Please pay attention to the sign of sigma in the output file.

    ymat_e=-(sigma_tot-2.0*vsigma_d0+2.0*vsigma_b0)/10.0-at.T
    if globalv.bulk_deform=='yes':
      #Please pay attention to the sign of sigma in the output file.
      # for ndm sigma_elastic=-sigma_output
      ymat=-((sigma_tot-sigma_tot_b)-2.0*vsigma_d0+2.0*vsigma_b0)/10.0

    Wfact = np.array([1., 100., 1., 1., 1., 1., 1.])
    weight = np.zeros(len(ymat))
    for i in range(len(ymat)):
      weight[i] = Wfact[i/globalv.nt_strain]
      # remove middle points for every strain type:
      # if ((i%nlen) > 3) and ((i%nlen) < 7):
      #  #weight[i] = 0.
      #  weight[i] = 1.0
    #debug print weight
    Wmat=np.diag(weight)
      #the same thing using my formula
    wp_e=np.matmul(np.matmul(np.linalg.inv(np.matmul(strain_tot.T,np.matmul(Wmat,strain_tot))),np.matmul(strain_tot.T,Wmat)),ymat_e)
    print("maison: alpha_ijkl using C_ijkl: polarizibility matrix in Voigt notation (eV):")
    alpha0=-wp_e*V0b/globalv.unit_conversion1
    printMatrixE(alpha0)


    if globalv.bulk_deform=='yes':

      wp=np.matmul(np.matmul(np.linalg.inv(np.matmul(strain_tot.T,np.matmul(Wmat,strain_tot))),np.matmul(strain_tot.T,Wmat)),ymat)


      print("maison: alpha_ijkl without C_ijkl: polarizibility matrix in Voigt notation (eV):")
      alpha = (-wp) *V0b/globalv.unit_conversion1
      printMatrixE(alpha)


    if globalv.bulk_deform=='yes':
      return  alpha
    else:
      return  alpha0




def fit_sigma_lmfit (sigma_tot, sigma_tot_b, strain_tot, vsigma_d0, vsigma_b0, V0d, V0b, cij):

    at=np.matmul(cij.T,strain_tot.T)

    #Please pay attention to the sign of sigma in the output file.

    ymat_e=-(sigma_tot-2.0*vsigma_d0+2.0*vsigma_b0)/10.0-at.T
    if globalv.bulk_deform=='yes':
      #Please pay attention to the sign of sigma in the output file.
      # for ndm sigma_elastic=-sigma_output
      ymat=-((sigma_tot-sigma_tot_b)-2.0*vsigma_d0+2.0*vsigma_b0)/10.0








# Denise --------------------------------------------------------------



def fit_scipy_leastq (ene_tot, ene_tot_b, ene_tot_corr, ecorr_sym, strain_tot, Pij0, nat0b, nat0d, E0b, E0d, V0b, V0d, cij, alpha):


  unit_conversion1 = globalv.unit_conversion1

  al=np.matmul(strain_tot,np.matmul(cij.T,strain_tot.T))*V0b/(2.0*unit_conversion1)

  # the elastic bulk energy is accounted only from elastic constants ... (77,1) object
  y_atomic=ene_tot[:,0]-E0d-np.diagonal(al)

  imiddle=int((globalv.nt_strain-1)/2)
  if globalv.bulk_deform=='yes':
     #the elastic energy of bulk deformation is computed atomisticaly ene_i_b - E0b
     y_atomic_full=ene_tot[:,0]-E0d \
                   -float(nat0d)/float(nat0b)*(ene_tot_b[:,0]-E0b)+ecorr_sym
     #the elastic energy of bulk deformation is computed atomisticaly ene_i_b - E0b + aneto correction
     print(imiddle)
     y_atomic_full_corr=ene_tot[:,0] -E0d-ene_tot_corr[:,0]+ene_tot_corr[imiddle,0] \
                   -float(nat0d)/float(nat0b)*(ene_tot_b[:,0]-E0b)+ecorr_sym

     y_atomic_full_corr[imiddle]=y_atomic_full_corr[imiddle]-ecorr_sym
     y_atomic_full[imiddle]=y_atomic_full[imiddle]-ecorr_sym

  al0=np.matmul(-Pij0,strain_tot.T)
  y_data=y_atomic_full_corr-al0[0,:]
  #debug print al0.shape, y_atomic_full_corr.shape, y_data.shape
  #al0=al
  #al0.shape 1,1 , Pij0.shape 1,6, strain_i.shape 11,6
  aint=np.matmul(strain_tot,np.matmul(-alpha,strain_tot.T))/2.0
  #al1=np.matmul(-Pij0,strain_tot.T)-np.diagonal(aint)/2.0
  al1=np.diagonal(aint)
  #debug 77, print al1.shape

  #convert diagonal form of  2D array alpha(6,6)  into beta(21) 1D array
#old  ibeta, beta, ij_beta = pack_alpha_into_beta(-alpha)
#old  y_s=model(beta, strain_tot, ibeta)
#old  res=least_squares(fun,beta,jac=jac, args=(strain_tot, y_data, ibeta), verbose=1, gtol=1.e-9)


  #uncomment for SIA
  #ibeta, beta, ij_beta = pack_alpha_into_beta_min_sia(-alpha)
  #ibeta, beta, ij_beta = pack_alpha_into_beta_min_sia_x(-alpha) # TJ
  ibeta, beta, ij_beta =pack_alpha_into_beta_sad_vac(-alpha) # TJ: works also for SIA at saddle position
  y_s=model_new(beta, strain_tot, ij_beta)
  res=least_squares(fun_new,beta,jac=jac_new, args=(strain_tot, y_data, ij_beta), verbose=1, gtol=1.e-6)

  y_test=model_new(res.x,strain_tot,ij_beta)
  try:
    import pylab
    #pylab.plot(y_s,y_data,'o')
    pylab.plot(y_data,al1,'x')
    pylab.plot(y_data,y_data,'-')
    pylab.plot(y_data,y_test,'o')
    pylab.show()
  except:
    pass
  #old alpha_fin=pack_beta_into_alpha(res.x,ibeta)
  alpha_fin=pack_beta_into_alpha_new(res.x,ij_beta)
  print("scipy leatsq: alpha_ijkl without C_ijkl: polarizibility matrix in Voigt notation (eV):")
  printMatrixE(-alpha_fin)



  return -alpha_fin

def pack_alpha_into_beta (alpha):

  beta = np.zeros(21)
  ibeta=np.array([]).reshape(0,2)
  ij_beta=[]

  kk = 0
  for ii in range(6):
    for jj in range(6):
      if jj > ii :
        beta[kk] = alpha[ii,jj]
        ij_beta_temp=[[ii,jj],[jj,ii]]
        ij_beta.append(ij_beta_temp)
        kk = kk+1
        ibeta=np.vstack([ibeta,np.array([ii,jj])])
      elif jj==ii :
        ij_beta_tmp=[[ii,jj]]
        ij_beta.append(ij_beta_tmp)
        beta[kk] = alpha[ii,jj]
        kk = kk+1
        ibeta=np.vstack([ibeta,np.array([ii,jj])])

  return ibeta.astype(int), beta, ij_beta


def pack_alpha_into_beta_sad_vac (alpha): # TJ: jump from [100] to [010] (saddle point direction = [110]); works also for SIA

  beta = np.zeros(9)
  ij_beta=[]
  #start with the diagonal ...
  beta[0]=alpha[0,0]
  ij_tmp=[[0,0],[1,1]]
  ij_beta.append(ij_tmp)
  beta[1]=alpha[2,2]
  ij_tmp=[[2,2]]
  ij_beta.append(ij_tmp)

  beta[2]=alpha[3,3]
  ij_tmp=[[3,3],[4,4]]
  ij_beta.append(ij_tmp)
  beta[3]=alpha[5,5]
  ij_tmp=[[5,5]]
  ij_beta.append(ij_tmp)



  beta[4]=alpha[0,1]
  ij_tmp=[[0,1],[1,0]]
  ij_beta.append(ij_tmp)
  beta[5]=alpha[1,2]
  ij_tmp=[[0,2],[1,2],[2,0],[2,1]]
  ij_beta.append(ij_tmp)


  beta[6]=alpha[0,5]
  ij_tmp=[[0,5],[1,5],[5,0],[5,1]]
  ij_beta.append(ij_tmp)
  beta[7]=alpha[2,5]
  ij_tmp=[[2,5],[5,2]]
  ij_beta.append(ij_tmp)

  beta[8]=alpha[3,4]
  ij_tmp=[[3,4],[4,3]]
  ij_beta.append(ij_tmp)



  ibeta=np.array([]).reshape(0,2)

  return ibeta.astype(int), beta,ij_beta

def pack_alpha_into_beta_min_vac (alpha):

  beta = np.zeros(3)
  ij_beta=[]
  beta[0]=alpha[0,0]
  ij_tmp=[[0,0],[1,1],[2,2]]
  ij_beta.append(ij_tmp)
  beta[1]=alpha[1,2]
  ij_tmp=[[0,1],[0,2],[1,2],[1,0],[2,0],[2,1]]
  ij_beta.append(ij_tmp)
  beta[2]=alpha[3,3]
  ij_tmp=[[3,3],[4,4],[5,5]]
  ij_beta.append(ij_tmp)
  ibeta=np.array([]).reshape(0,2)

  return ibeta.astype(int), beta,ij_beta


def pack_alpha_into_beta_min_sia (alpha):

  beta = np.zeros(6)
  ij_beta=[]
  beta[0]=alpha[0,0]
  ij_tmp=[[0,0],[1,1]]
  ij_beta.append(ij_tmp)

  beta[1]=alpha[2,2]
  ij_tmp=[[2,2]]
  ij_beta.append(ij_tmp)

  beta[2]=alpha[0,1]
  ij_tmp=[[0,1],[1,0]]
  ij_beta.append(ij_tmp)


  beta[3]=alpha[0,2]
  ij_tmp=[[0,2],[1,2], [2,0], [2,1]]
  ij_beta.append(ij_tmp)

  beta[4]=alpha[3,3]
  ij_tmp=[[3,3], [4,4]]
  ij_beta.append(ij_tmp)

  beta[5]=alpha[5,5]
  ij_tmp=[[5,5]]
  ij_beta.append(ij_tmp)

  ibeta=np.array([]).reshape(0,2)

  return ibeta.astype(int), beta,ij_beta

def pack_alpha_into_beta_min_sia_x (alpha): # TJ: SIA oriented along [100] = e_x

  beta = np.zeros(6)
  ij_beta=[]
  beta[0]=alpha[0,0]
  ij_tmp=[[0,0]]
  ij_beta.append(ij_tmp)

  beta[1]=alpha[2,2]
  ij_tmp=[[1,1],[2,2]]
  ij_beta.append(ij_tmp)

  beta[2]=alpha[0,1]
  ij_tmp=[[0,1],[1,0],[2,0],[0,2]]
  ij_beta.append(ij_tmp)

  beta[3]=alpha[1,2]
  ij_tmp=[[1,2],[2,1]]
  ij_beta.append(ij_tmp)

  beta[4]=alpha[3,3]
  ij_tmp=[[3,3]]
  ij_beta.append(ij_tmp)

  beta[5]=alpha[5,5]
  ij_tmp=[[4,4],[5,5]]
  ij_beta.append(ij_tmp)

  ibeta=np.array([]).reshape(0,2)

  return ibeta.astype(int), beta,ij_beta


def pack_beta_into_alpha(beta, ibeta):

  alpha = np.zeros((6,6))
  for i in range(len(beta)):
    if ibeta[i,0]==ibeta[i,1]:
      alpha[ibeta[i,0],ibeta[i,1]]=beta[i]
    else:
      alpha[ibeta[i,0],ibeta[i,1]]=beta[i]/2.0
      alpha[ibeta[i,1],ibeta[i,0]]=beta[i]/2.0

  return alpha

def pack_beta_into_alpha_new(beta, ij_beta):

  alpha = np.zeros((6,6))
  for kk in range(len(beta)):
       for i in range(len(ij_beta[kk])):
            alpha[ij_beta[kk][i][0],ij_beta[kk][i][1]]=beta[kk]

  return alpha







def fun(beta,x,y,ibeta):

   return model(beta,x,ibeta)-y

def fun_new(beta,x,y,ij_beta):

   return model_new(beta,x,ij_beta)-y



def  model(beta,x,ibeta):
   yl=np.zeros(x[:,0].size)
   for i in range(len(beta)):
    yl += beta[i]*x[:,ibeta[i,0]]*x[:,ibeta[i,1]]
   return yl/2.0

def  model_new(beta,x,ij_beta):
   yl=np.zeros(x[:,0].size)
   for k in range(len(beta)):
      tmp=0
      for j in range(len(ij_beta[k])):
        tmp += x[:,ij_beta[k][j][0]]*x[:,ij_beta[k][j][1]]
      yl += beta[k]*tmp
   return yl/2.0

def jac(beta,x,y,ibeta):

   J=np.empty((y.size,beta.size))
   for k in range(beta.size):
     J[:,k] = x[:,ibeta[k,0]]*x[:,ibeta[k,1]]
   return J

def jac_new(beta,x,y,ij_beta):

   J=np.empty((y.size,beta.size))
   for k in range(beta.size):
      tmp=0
      for j in range(len(ij_beta[k])):
        tmp += x[:,ij_beta[k][j][0]]*x[:,ij_beta[k][j][1]]
      J[:,k] = tmp
   return J




def func_ecorr_sym (def_type):

  if def_type==2:
     #vac_sad
     ecorr_sym=0.00075
     #vac_min ecorr_sym=0.0
     #sia min corr ecorr_sym=0.00085

  if def_type==3:
     #vac_sad
     ecorr_sym=0.00075
     #vac_min ecorr_sym=0.0
     #sia min corr ecorr_sym=0.00085


  if def_type==4:
     #vac_sad
     ecorr_sym=0.0006
     #vac min ecorr_sym=0.0003
     #sia min corr
  if def_type==5:
     #vac_sad
     ecorr_sym=0.00054
     #ecorr_sym=0.0
     #vac_min ecorr_sym=0.0003
  if def_type==6:
     #vac_sad
     ecorr_sym=0.0006
     #ecorr_sym=0.0
     #vac_min ecorr_sym=0.0003
  if def_type==7:
     #vac_sad
     ecorr_sym=0.0006
     #ecorr_sym=0.0
     #vac_min ecorr_sym=0.0003

  if def_type==8:
     #vac_sad
     ecorr_sym=0.0006
     #ecorr_sym=0.0
     #vac_min ecorr_sym=0.0003
  if def_type==9:
     #vac_sad ecorr_sym=0.0006
     ecorr_sym=0.0
     #vac_min ecorr_sym=0.0003

  if def_type==23:
     #vac_sad ecorr_sym=0.0006
     ecorr_sym=0.0006
     #vac_min ecorr_sym=0.0003
  ecorr_sym=0.0


  return ecorr_sym
