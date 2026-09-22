import globalv
import numpy as np

def voigt_elastic_matrix_bcc(c11,c12,c44):

   voigt=np.matrix(np.zeros((6,6)))
   
   voigt[0,0]=float(c11)
   voigt[0,1]=float(c12)
   voigt[1,0]=float(c12)
   voigt[0,2]=float(c12)
   voigt[2,0]=float(c12)
   voigt[1,1]=float(c11)
   voigt[1,2]=float(c12)
   voigt[2,1]=float(c12)
   voigt[2,2]=float(c11)
   voigt[3,3]=float(c44)
   voigt[4,4]=float(c44)
   voigt[5,5]=float(c44)

   return voigt;

def index_voigt(ii,jj):
   i=int(ii)
   j=int(jj)
   if (i>2) or (j>2):
      print('ERROR in index_voigt',i,j)
      exit(0)

   if i==j:
      ivoigt=i
   else:

      if ((i==1 and j==2) or (i==2 and j==1)): 
        ivoigt=3 
      if ((i==0 and j==2) or (i==2 and j==0)): 
        ivoigt=4
      if ((i==0 and j==1) or (i==1 and j==0)):
        ivoigt=5
   return int(ivoigt);

def factor_voigt(ii):
   i=int(ii)

   if (i==1 or i==2 or i==0): 
     factor=1.0 
   if (i==3 or i==4 or i==5): 
    factor=2.0 
   return float(factor);



def deform_matrix(strain_tensor):
  e1,e2,e3,e4,e5,e6 = strain_tensor
#  e_p = np.mat ([[(1.0+e1),0.5*e6,0.5*e5],\
#                  [0.5*e6,(1.0+e2),0.5*e4],\
#                  [0.5*e5,0.5*e4,(1.0+e3)]])

  e_p = np.mat ([[(0.0+e1),0.5*e6,0.5*e5],\
                  [0.5*e6,(0.0+e2),0.5*e4],\
                  [0.5*e5,0.5*e4,(0.0+e3)]])



  return e_p;


def pack_sigma_in_vsigma(sigma):
  vsigma=np.array([[sigma[0,0],sigma[1,1],sigma[2,2],\
                    0.5*(sigma[1,2]+sigma[2,1]),\
                    0.5*(sigma[0,2]+sigma[2,0]),\
                    0.5*(sigma[1,0]+sigma[0,1])]])
  return vsigma;
def pack_strain_in_vstrain(sigma):
  vsigma=np.array([[sigma[0,0],sigma[1,1],sigma[2,2],\
                    1.0*(sigma[1,2]+sigma[2,1]),\
                    1.0*(sigma[0,2]+sigma[2,0]),\
                    1.0*(sigma[1,0]+sigma[0,1])]])
  return vsigma;

def print_strain(t_ini,t_fin,npoints, i,t,strain_tensor):
   fs=open('strain_save.dat','w')

   fs.write('%s \n'%t_ini)
   fs.write('%s \n'%t_fin)
   fs.write('%s \n'%npoints)
   fs.write('%s \n'%i)
   fs.write('%s \n'%t)
   for i in range(len(strain_tensor)):
      fs.write('%20.15E \n'%strain_tensor[i])
   fs.close
   return;


