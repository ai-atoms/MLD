import re
import os
import math
import matplotlib
from matplotlib import pyplot as pl
import numpy as np
from scipy.interpolate import CubicSpline

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
   print("get_sigma_from_lammps: error in sigma. Probably, the output of LAMMPS is not patched")
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
   print('get_energy_from_ndm: There is no energy output in the directory %s'%os.getcwd())
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
     print("get_energy_from_lammps: error in energy. Probably, the output of LAMMPS is not patched")
     exit(0)
   #if math.fabs(energy - p_energy) > 1.e-5 :
   # print 'get_energy_from_lammps: incorrect or inconsistency in energy'
   # exit(0)

 return p_energy;


def get_volume_from_lammps(outNDM, Patched=None):
 f=open(outNDM,'r')
 word=[]

 if Patched is None:
   print('get_volume_from_lammps: The output %s should pe patched with Volume keyword'%os.getcwd())
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
     print("get_volume_from_lammps: error in volume. Probably, the output of LAMMPS is not patched")
     exit(0)

 return float(volume);



def get_climb_replica_from_lammps(outNDM, Patched=None):
 f=open(outNDM,'r')
 word=[]

 if Patched is None:
   print('get_climb_replica_from_lammps: The output %s should comes from NEB'%os.getcwd())
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
     print("get_climb_replica_from_lammps: Probably, the output of LAMMPS not comes from NEB")
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
  ax.plot(xepsilon,y_atomic_smooth(xepsilon),linestyle='--',markersize=6,color='slateblue', label='interpolation')

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
line=return_last_line('out.run')
word=line.split()
line_last=np.asarray(line.split()[9:], dtype=float).reshape(-1,2)
print(line_last)
rcoord=line_last[:,0]
energy=line_last[:,1]
print(rcoord)
print(energy)


energy_new=[]
for i in range(len(energy)):
      ene_tmp=get_energy_from_lammps('screen.%i'%i, Patched='Yes')
      energy_new.append(float(ene_tmp))



f=open('ene_path','w')

for i in range(len(energy)):
      print(rcoord[i], energy[i], energy_new[i])
      f.write('%s %s \%s\n'%(rcoord[i], energy_new[i]-energy_new[0],    energy_new[i]))



xp_vol=np.array(rcoord).astype(float)
y_atomic=np.array(energy_new).astype(float)

plot_data (xp_vol, y_atomic)

ymin=np.amin(y_atomic)
min_loc=np.argmin(y_atomic)

ymax=np.amax(y_atomic)
max_loc=np.argmax(y_atomic)

if (min_loc != 0):
  print ("NEB WARNING the 0 image is not the lowest energy along the path")

isaddle=get_climb_replica_from_lammps('out.run', Patched='Yes')

file_saddle='screen.%i'%isaddle
print(file_saddle)

ene_saddle=get_energy_from_lammps('screen.%i'%isaddle, Patched='Yes')
ene_min=get_energy_from_lammps('screen.0', Patched='Yes')
sigma=get_sigma_from_lammps('screen.%i'%isaddle, Patched='Yes')



print('From out_neb', ymin, ymax, ymax-ymin)
print('From out_neb', min_loc, max_loc)
print('From energy ', ene_min, ene_saddle, ene_saddle-ene_min)
print('From energy ', isaddle)


print(sigma)
