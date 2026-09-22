import array
import math
import os 
import sys
import numpy as np
from tempfile  import mkstemp
from shutil    import move





def read_eigenvalues_from_ndm(file_eigenvalues):
  icount=0
  omega=[]
  feigen=open(file_eigenvalues,'r')
  lines_feigen=feigen.read().splitlines()
  for i in  range(len(lines_feigen)):
     if not lines_feigen[i]=='':
       word=lines_feigen[i].split()
       if not (word[0]=='#'):
         icount=icount+1
         omega.append(float(word[1]))
  omega=np.array(omega)
#  print icount , omega.size
  return omega

def read_eigenvalues_from_vasp(file_eigenvalues):
  icount=0
  omega=[]
  feigen=open(file_eigenvalues,'r')
  lines_feigen=feigen.read().splitlines()
  for i in  range(len(lines_feigen)):
     if not lines_feigen[i]=='':
       word=lines_feigen[i].split()
       if not (word[0]=='#'):
         icount=icount+1
         omega.append(float(word[0]))
  omega=np.array(omega)
#debug   print icount , omega.size
  return omega



def free_energy_classical(omega,limit_omega,temperature):
   thz_to_ev=0.004135665538536
   temperature_to_ev=8.6173303E-05 
   tev=temperature*temperature_to_ev
   xp=omega[:]*thz_to_ev/tev

   free_energy = np.sum(tev*np.log(xp))
   entropy = np.sum(temperature_to_ev*(1.0 - np.log(xp)))
   cv= temperature_to_ev * len(omega[:])
   internal_energy = tev * len(omega[:])

   return free_energy, entropy, cv, internal_energy




def free_energy_quantum(omega,limit_omega,temperature):
   thz_to_ev=0.004135665538536
   temperature_to_ev=8.6173303E-05 

   tev=temperature*temperature_to_ev
   xp2 = omega[:]*thz_to_ev/(2.0*tev)

   free_energy  = tev*np.sum(np.log(2.0*np.sinh(xp2)))
   entropy  = temperature_to_ev*np.sum(-np.log(2.0*np.sinh(xp2))+xp2/np.tanh(xp2))
   cv = temperature_to_ev*np.sum((xp2**2/np.tanh(xp2)**2-xp2**2))
   internal_energy  = np.sum(omega/(np.tanh(xp2)*2.0))
   return free_energy, entropy, cv, internal_energy




def free_energy_quantum_old_slow(omega,limit_omega,temperature):
   thz_to_ev=0.004135665538536
   temperature_to_ev=8.6173303E-05 
   free_energy=0.0
   entropy=0.0
   cv=0.0
   internal_energy=0.0
   for i in  range(omega.size):
      fo=0.0
      so=0.0
      dso=0.0
      uo=0.0
      if np.fabs(omega[i]) >= limit_omega:
         if (omega[i]<0.0):
           print('WARNING imaginary freq...in free energy quantum', i, omega[i])
           os.system('pwd')
         omv=omega[i]*thz_to_ev
         tev=temperature*temperature_to_ev
         xp2=omv/(tev*2.0)
         xsinh = np.sinh(xp2)
         xtanh = np.tanh(xp2)
         tempx=1.0-np.exp(-omv/tev)
         #tempy=math.exp(omv/tev)-1.0
         #free energy F
         #old version fo = omv/2.0 + tev*math.log(tempx)
         fo = tev*np.log(2.0*xsinh)
         # S entropy (ev/K units) ... so TS with T in K is eV 
         #old version so = temperature_to_ev*(-math.log(tempx)+omv/(tev*tempy))
         so = temperature_to_ev*(-np.log(2.0*xsinh)+xp2/xtanh)
         dso =temperature_to_ev*(xp2**2/xtanh**2-xp2**2)
         # internal energy ...
         #old version uo = omv/2.0 + omv/tempy
         uo=omv/(xtanh*2.0)

      free_energy=free_energy+fo
      entropy=entropy+so
      cv=cv+dso
      internal_energy=internal_energy+uo
   return free_energy, entropy, cv, internal_energy

def free_energy_classical_old_slow(omega,limit_omega,temperature):
   thz_to_ev=0.004135665538536
   temperature_to_ev=8.6173303E-05 
   free_energy=0.0
   entropy=0.0
   cv=0.0
   internal_energy=0.0
   icount=0
   for i in  range(omega.size):
      fo=0.0
      so=0.0
      dso=0.0
      uo=0.0
      if math.fabs(omega[i]) >= limit_omega:
         if (omega[i]<0.0):
           print('WARNING imaginary freq...in free_energy_classical', i,omega [i])
           os.system('pwd')
         omv=omega[i]*thz_to_ev
         tev=temperature*temperature_to_ev
         xp=omv/tev
         #free energy F
         fo = tev*math.log(xp)
         # S entropy (ev/K units) ... so TS xith T in K is eV 
         so = temperature_to_ev*(1.0 - math.log(xp))
         dso= temperature_to_ev
         # internal energy ...
         uo = tev
         icount=icount+1
      free_energy=free_energy+fo
      entropy=entropy+so
      cv=cv+dso
      internal_energy=internal_energy+uo

   return free_energy, entropy, cv, internal_energy






def get_roots_of_the_polynom_1D(val_min,val_max,poly_data, temperature):
   #val_max=val_max-0.3   
   pol_roots=np.roots(poly_data)
   itest=0
   #debug print val_min, val_max, pol_roots
   sol=[]
   for i in range(len(poly_data)-1):
      if abs(pol_roots[i].imag) < 1e-16:
        if (pol_roots[i] >= val_min) and (pol_roots[i] <= val_max):
          find_root=pol_roots[i].real
          sol.append(pol_roots[i])
          itest=itest+1

   if (itest==0):
      print('get_roots_of_the_poly: Fatal error there is no solution for this equation')
      print(val_min, val_max, pol_roots)
      print('temperature %s'%temperature)
      exit(0)
   if (itest>1):
      print('WARNING:  maybe more than one solution for the polynomial equation %s'%(itest))
      print(val_min, val_max, pol_roots)
      print('temperature %s'%temperature)
   return find_root,sol


def formation_energy(nbulk,ntotal,ebulk,etotal):
   ene = etotal - float(ntotal)/float(nbulk)*ebulk
   return ene;

def is_upper_triangular(mat):
    """test if 3x3 matrix is upper triangular"""

    def near0(x):
        """Test if a float is within .00001 of 0"""
        return abs(x) < 0.00001

    return near0(mat[1, 0]) and near0(mat[2, 0]) and near0(mat[2, 1]);


def right_hand_basis(A,B,C):
   """test if the reper is right handed"""

   return np.dot(np.cross(A,B),C) > 0 ; 

def tail(f, window):
    """ 
    Returns the last `window` lines of file `f` as a list.
    """
    if window == 0:
        return []
    BUFSIZ = 1024
    f.seek(0, 2)
    bytes = f.tell()
    size = window + 1 
    block = -1
    data = []
    while size > 0 and bytes > 0:
        if bytes - BUFSIZ > 0:
            # Seek back one whole BUFSIZ
            f.seek(block * BUFSIZ, 2)
            # read BUFFER
            data.insert(0, f.read(BUFSIZ))
        else:
            # file too small, start from begining
            f.seek(0,0)
            # only read what was not read
            data.insert(0, f.read(bytes))
        linesFound = data[0].count('\n')
        size -= linesFound
        bytes -= BUFSIZ
        block -= 1
    return ''.join(data).splitlines()[-window:]


def head(fname,nline):

   f=open(fname,'r')

   head=[next(f) for x in range(nline) ]
   f.close()

   return ''.join(data).splitlines()


def headn(file_name, n): 
    result = []
    nlines = 0 
    assert n >= 1
    for line in open(file_name):
        result.append(line)
        nlines += 1
        if nlines >= n:
            break
    return result




def printMatrixE(a):
   print("Matrix["+("%d" %a.shape[0])+"]["+("%d" %a.shape[1])+"]")
   rows = a.shape[0]
   cols = a.shape[1]
   for i in range(0,rows):
      for j in range(0,cols):
         print(("%8.3f" %a[i,j]), end=' ')
      print()
   print()

  

def search_word_replace_line(file_path, pattern, subst):
#from tempfile  import mkstemp
#from shutil    import move

   file_path = os.path.abspath(file_path)
   #Create temp file
   fh, abs_path = mkstemp()
   new_file = open(abs_path,'w')
   old_file = open(file_path)
   for line in old_file:
     if pattern in line:
        new_file.write(subst)
     else:
        new_file.write(line)
   #close temp file
   new_file.close()
   os.close(fh)
   old_file.close()
   #Remove original file
   os.remove(file_path)
   #Move new file
   move(abs_path, file_path)

