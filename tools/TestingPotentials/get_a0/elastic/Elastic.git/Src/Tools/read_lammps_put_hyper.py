import re
import os
import math
import matplotlib
from matplotlib import pyplot as pl
import numpy as np
from scipy.interpolate import CubicSpline



def read_lammps_data(dirSTR,name_of_input_file=None):

 if name_of_input_file is None:
   sgin= dirSTR+'/'+ 'structure.lmp'
 else: 
   sgin= dirSTR+'/'+ name_of_input_file
  
 if not os.path.exists(sgin):
    print "lammps file %s doesn't exist. Put the correct path"%str(sgin)
    exit(0)
 
 fgin=open(sgin,'r')
 lines_fgin=fgin.read().splitlines()
 icount=0
 itest=99999999999999999999999
 
 for i in range(len(lines_fgin)):
   if not lines_fgin[i]=='':
     word=lines_fgin[i].split()
     if i==2:
       if word[-1]=='atoms':
          nat_per_box=int(word[0])
          rcell=np.matrix(np.zeros((nat_per_box,3)))
          cell=np.matrix(np.zeros((nat_per_box,3)))
          itype=np.array(np.zeros(nat_per_box))
       else:
          print 'error in line 3 of LAMMPS'
          exit(0)

     if i==3:
       if word[-1]=='types':
         ntype=int(word[0])
       else:
         print 'error in line 4 of LAMMPS'
         exit(0)


     if i==5:
        xlo=float(word[0])
        xhi=float(word[1])
        if not(word[-1]=='xhi'):
          print 'not legal format for xhi'
          exit(0)

     if i==6:
        ylo=float(word[0])
        yhi=float(word[1])
        if not(word[-1]=='yhi'):
          print 'not legal format for yhi'
          exit(0)

     if i==7:
        zlo=float(word[0])
        zhi=float(word[1])
        if not(word[-1]=='zhi'):
          print 'not legal format for zhi'
          exit(0)

     if (i==8):
        xy=float(word[0])
        xz=float(word[1])
        yz=float(word[2])
        if not(word[-1]=='yz'):
          print 'not legal format for yz'
          exit(0)
        if (xlo > 0.0) or (ylo > 0.0) or (zlo > 0.0):
           print 'format unknown for xlo , ylo or zlo > 0'
           exit(0)
        #lammps_cell[0, 0] = xhi 
        #lammps_cell[1, 1] = yhi 
        #lammps_cell[2, 2] = zhi 
        #lammps_cell[0, 1] = xy  
        #lammps_cell[0, 2] = xz  
        #lammps_cell[1, 2] = yz  
        a,b,c,alpha,beta,gamma=convert_xyzlo_into_abc(xlo,ylo,zlo,xhi,yhi,zhi,xy,xz,yz)
        boxini= convert_abc_into_box_cell(a,b,c,alpha,beta,gamma)
     if i==10:
       if word[0]=='Masses':
          itest=16
       elif word[0]=='Atoms':
          itest=12
       else: 
          print 'error in LAMMPS file at line 10'
          exit(0)

     if (i>=itest) and (i <= (itest+nat_per_box)):
         iatom=int(word[0])
         itype[iatom-1]=int(word[1])
         rcell[iatom-1,0]=float(word[2])
         rcell[iatom-1,1]=float(word[3])
         rcell[iatom-1,2]=float(word[4])



 fgin.close()
 boxini_inv=np.linalg.inv(boxini)
 cell = np.matmul(boxini_inv,rcell.T).T

 return nat_per_box, boxini, itype, cell, rcell; 



def  convert_xyzlo_into_abc(xlo,ylo,zlo,xhi,yhi,zhi,xy,xz,yz):
# This function converts xlo,ylo... lammps into
# a, b, c, alpha, beta, gamma format.


 x_t = xhi-xlo
 y_t = yhi-ylo
 z_t = zhi-zlo

 a = x_t
 b = np.sqrt(y_t**2+xy**2)
 c = np.sqrt(z_t**2 + xz**2 + yz**2)
 alpha = np.arccos((xy*xz+y_t*yz)/(b*c))
 beta  = np.arccos(xz/c)
 gamma=np.arccos(xy/b)
 return a,b,c,alpha,beta,gamma;



def  convert_abc_into_box_cell(a_size,b_size,c_size,alpha,beta,gamma):

#  This function converts conventional LAMMPS vectors
#   defined by a b c alpha beta gamma,
#   into a lower triangular matrix:
#         |  H(0,0)    0       0     |
#    H =  |  H(1,0)  H(1,1)    0     |
#         |  H(2,0)  H(2,1)  H(2,2)  |
#The matrix which pre-multiplies the column vector of the 
#fractional crystallographic coordinates to yield the distributed 
#coordinates in the a, b, b system is:

#      a   b(cos(gamma))   c(cos(beta))
#      0   b(sin(gamma))   c(cos(alpha) - cos(beta) cos(gamma)) / sin(gamma)
#      0   0               V/(ab sin(gamma))

#V = abc(1 - cos**2(alpha) - cos**2(beta) - cos**2(gamma) + 2(cos(alpha) cos(beta) cos(gamma)))**1/2 

 box_cell=np.zeros((3,3))

 #a_x
 box_cell[0,0] = a_size
 #b_x
 box_cell[1,0] = b_size*np.cos(gamma)
 #b_y
 box_cell[1,1] = b_size*np.sin(gamma)
 #c_x
 box_cell[2,0] = c_size*np.cos(beta)
 #c_y
 box_cell[2,1] = c_size*( ( np.cos(alpha)-np.cos(beta)*np.cos(gamma) )/ \
                         np.sin(gamma) )
 #c_z
 box_cell[2,2] = c_size*np.sqrt(                                    \
                       np.sin(gamma)**2                              \
                      -np.cos(beta)**2 - np.cos(alpha)**2            \
                      +2.0*np.cos(alpha)*np.cos(beta)*np.cos(gamma) \
                     )/np.sin(gamma)                  
                     
 return box_cell; 





def is_upper_triangular(mat):
    """test if 3x3 matrix is upper triangular"""

    def near0(x):
        """Test if a float is within .00001 of 0"""
        return abs(x) < 0.00001

    return near0(mat[1, 0]) and near0(mat[2, 0]) and near0(mat[2, 1])


def right_hand_basis(A,B,C):
   """test if the reper is right handed"""

   return np.dot(np.cross(A,B),C) > 0
   

def convert_cell_from_gin_to_lammps(box_cell):
    """
    Convert a parallelpiped (forming right hand basis)
    to lower triangular matrix LAMMPS can accept. This
    function transposes cell matrix so the bases are column vectors
    """
    cell = np.matrix.transpose(box_cell)

    if not is_upper_triangular(cell):
        # rotate bases into triangular matrix
        tri_mat = np.zeros((3, 3))
        A = np.array(cell[:, 0])
        B = np.array(cell[:, 1])
        C = np.array(cell[:, 2])

        #print 'A', A, np.shape(A)
        #print 'B', B, np.shape(B)
        #print 'C', C, np.shape(C)
        if not right_hand_basis(A,B,C):
             print"WARNING: your reper is not right handed."
             print"WARNING: This is a critical issue. The LAMMPS results are wrong !!!!!" 


        tri_mat[0, 0] = norm(A)
        Ahat = A / norm(A)
        AxBhat = np.cross(A, B) / norm(np.cross(A, B))
        tri_mat[0, 1] = np.dot(B, Ahat)
        tri_mat[1, 1] = norm(np.cross(Ahat, B))
        tri_mat[0, 2] = np.dot(C, Ahat)
        tri_mat[1, 2] = np.dot(C, np.cross(AxBhat, Ahat))
        tri_mat[2, 2] = norm(np.dot(C, AxBhat))

        # create and save the transformation for coordinates
        volume = np.linalg.det(box_cell)
        trans = np.array([np.cross(B, C), np.cross(C, A), np.cross(A, B)])
        trans = trans / volume
        coord_transform = np.dot(tri_mat , trans)

        return tri_mat, coord_transform
    else:
        return cell, np.eye(3)


#interaction with LAMMPS .....







#LAMMPS related ..............











def get_sigma_from_lammps(outNDM, Patched=None):
 """ 
 Get sigma from patched output of LAMMPS
 patched means that the component of the stress are labeled by sig_xx, sig_xy, sig_xz, sig_yz, sig_zz
 output is in kbar 
 """

 f=open(outNDM,'r')
 word=[]
 i=0
 sigma=np.zeros((3,3))
 i=0
 for line in f:
   if re.match('(.*)sig_xx(.*)', line):
     word=line.split()
     sigma[0,0]=float(word[-1])
     i+=1
   if re.match('(.*)sig_yy(.*)', line):
     word=line.split()
     sigma[1,1]=float(word[-1])
     i+=1
   if re.match('(.*)sig_xy(.*)', line):
     word=line.split()
     sigma[0,1]=float(word[-1])
     sigma[1,0]=float(word[-1])
     i+=1
   if re.match('(.*)sig_xz(.*)', line):
     word=line.split()
     sigma[0,2]=float(word[-1])
     sigma[2,0]=float(word[-1])
     i+=1
   if re.match('(.*)sig_yz(.*)', line):
     word=line.split()
     sigma[1,2]=float(word[-1])
     sigma[2,1]=float(word[-1])
     i+=1
   if re.match('(.*)sig_zz(.*)', line):
     word=line.split()
     sigma[2,2]=float(word[-1])
     i+=1


 if i!=6:
   print"get_sigma_from_lammps: error in sigma. Probably, the output of LAMMPS is not patched"
   exit(0)
  
#patched version ...convert in kbar:
 sigma=sigma*1.e-3 
 return sigma;

def get_energy_from_lammps(outNDM, Patched=None):
 f=open(outNDM,'r')
 word=[]
 i=0
 for line in f:
   if i==1:
     word=line.split()
     energy=float(word[-1])
     break
   if re.match('(.*) Energy initial(.*)', line):
     word=line.split()
     i=i+1

 if len(word)==0:
   print 'get_energy_from_ndm: There is no energy output in the directory %s'%os.getcwd()
   exit(0)


#patched version ...

 if Patched is not None: 
   f.seek(0)
   for line in f:
     if re.match('(.*)Energy_box(.*)', line):
       word=line.split()
       p_energy=float(word[-1])
       break

   f.close()
   
   try:
       p_energy
   except NameError: 
     print"get_energy_from_lammps: error in energy. Probably, the output of LAMMPS is not patched"
     exit(0)
   #if math.fabs(energy - p_energy) > 1.e-5 :
   # print 'get_energy_from_lammps: incorrect or inconsistency in energy'
   # exit(0) 

 return p_energy;


def get_volume_from_lammps(outNDM, Patched=None):
 f=open(outNDM,'r')
 word=[]
 
 if Patched is None: 
   print 'get_volume_from_lammps: The output %s should pe patched with Volume keyword'%os.getcwd()
   exit(0)

#patched version ...

 if Patched is not None: 
   f.seek(0)
   for line in f:
     if re.match('(.*)Volume(.*)', line):
       word=line.split()
       volume=float(word[-1])
       break

   f.close()
   
   try:
       volume
   except NameError: 
     print"get_volume_from_lammps: error in volume. Probably, the output of LAMMPS is not patched"
     exit(0)
   
 return float(volume);



def get_climb_replica_from_lammps(outNDM, Patched=None):
 f=open(outNDM,'r')
 word=[]
 
 if Patched is None: 
   print 'get_climb_replica_from_lammps: The output %s should comes from NEB'%os.getcwd()
   exit(0)

#patched version ...

 if Patched is not None: 
   f.seek(0)
   for line in f:
     if re.match('(.*)Climbing replica(.*)', line):
       word=line.split()
       no_replica=float(word[-1])
       break

   f.close()
   
   try:
       no_replica
   except NameError: 
     print"get_climb_replica_from_lammps: Probably, the output of LAMMPS not comes from NEB"
     exit(0)
   
 return float(no_replica);


def plot_data(xp_vol, y_atomic):

  ymin=np.amin(y_atomic)
  y_tmp=y_atomic-ymin
  y_atomic=y_tmp
  xmin=xp_vol[0]-0.05/100.0
  xmax=xp_vol[-1]+0.05/100.0
  xepsilon=np.linspace(xmin,xmax,num=200,endpoint=True)
  #y_atomic_smooth=interp1d(xp_vol/100.0,y_atomic,kind='cubic')
  y_atomic_smooth=CubicSpline(xp_vol,y_atomic)
  fig=pl.figure()
  ax=fig.add_subplot(111)
  ax.plot(xp_vol,y_atomic,'o',markersize=6,color='slateblue',label='atomic') 
  ax.plot(xepsilon,y_atomic_smooth(xepsilon),linestyle='-',markersize=6,color='slateblue', label='interpolation') 

  ax.set_title('NEB energy reference %s eV'%ymin)
  ax.set_ylabel(r" $E (r)$ (eV) ", color='slateblue',fontsize='20')
  ax.set_xlabel(r" r (reaction coordinate) ")
  ax.set_xlim([xmin,xmax])
  ax.tick_params(axis='x',labelsize='16')
  ax.tick_params(axis='y',labelsize='16')
  ax.legend(loc='best') 
  defsize=fig.get_size_inches()
  fig.set_size_inches( (defsize[0]*1.5,defsize[1]*1.5) )
  fig.savefig('%s.png'%'graph')
  fig.savefig('%s.eps'%'graph')

  pl.show()

  return;



def return_last_line (outNDM):


  with open(outNDM,'r') as fneb:
    lines=fneb.read().splitlines()
    last_line=lines[-1]
  
  return last_line;


rcoord=[]
energy=[]
line=return_last_line('out_neb')
word=line.split()
itest=0
ilast=0
for i in range(len(word)):
   if float(word[i])==0:
     itest=1
   if (itest >= 1):
     if (itest-1)%2 == 1:
         energy.append(word[i])
         if ilast==1:
           break
     if (itest-1)%2 == 0:
         rcoord.append(word[i])
         if float(word[i]) == 1:
           ilast=1
     itest=itest+1


for i in range(len(energy)):
   print rcoord[i], energy[i]


xp_vol=np.array(rcoord).astype(np.float)
y_atomic=np.array(energy).astype(np.float)


ymin=np.amin(y_atomic)
min_loc=np.argmin(y_atomic)

ymax=np.amax(y_atomic)
max_loc=np.argmax(y_atomic)

if (min_loc != 0):
  print ("NEB WARNING the 0 image is not the lowest energy along the path") 
 
isaddle=get_climb_replica_from_lammps('out_neb', Patched='Yes')

file_saddle='screen.%i'%isaddle
print file_saddle

ene_saddle=get_energy_from_lammps('screen.%i'%isaddle, Patched='Yes')
ene_min=get_energy_from_lammps('screen.0', Patched='Yes')
sigma=get_sigma_from_lammps('screen.%i'%isaddle, Patched='Yes')



print 'From out_neb', ymin, ymax, ymax-ymin
print 'From out_neb', min_loc, max_loc
print 'From energy ', ene_min, ene_saddle, ene_saddle-ene_min
print 'From energy ', isaddle

nat_per_box, boxini, itype, cell, rcell=read_lammps_data('.','final_neb.0')
tri_mat, beta_matrix = convert_cell_from_gin_to_lammps(boxini)

xhi = tri_mat[0, 0]
yhi = tri_mat[1, 1]
zhi = tri_mat[2, 2]
xy = tri_mat[0, 1]
xz = tri_mat[0, 2]
yz = tri_mat[1, 2]
fi=open('SC_INFO','w')
fi.write('# images from %i neb images\n'%len(energy))
fi.write('0.0  %.14f\n'%xhi)
fi.write('0.0  %.14f\n'%yhi)
fi.write('0.0  %.14f\n'%zhi)
fi.write('  %.14f %.14f %.14f\n'%(xy, xz,yz))
fi.write('0.0  0.0\n')
for i in range(len(energy)):
      
      fi.write('_knot/PREFIX/knot_%i.xyz\n'%i)

fi.close()





for i in range(len(energy)):
   file_neb='final_neb.'+str(i)
   nat_per_box, boxini, itype, cell, rcell=read_lammps_data('.',file_neb)
   fi=open('knot_%i.xyz'%i,'w')
   for ia in range(nat_per_box):
      fi.write('%.15e %.15e %.15e\n'%(rcell[ia,0], rcell[ia,1], rcell[ia,2]))
   fi.close()





plot_data (xp_vol, y_atomic)

print sigma










