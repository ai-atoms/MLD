#!/usr/bin/python
import numpy as np
import os
from numpy.linalg import norm
import sys

def read_gin (dirSTR, ginFile=None):
#   root_dir=globalv.root_dir
#   dirSTR=root_dir+'/Structure'
     if ginFile is None:
       sgin= dirSTR+'/'+ 'structure.gin'
     else: 
       sgin= dirSTR+'/'+ ginFile
  
     if not os.path.exists(sgin):
        print("Gin file %s doesn't exist. Put the correct path in setup_ndm and structure"%str(sgin))
        exit(0)
  
     # should be readed in gin ...
     fgin=open(sgin,'r')
     lines_fgin=fgin.read().splitlines()
     icount=0
     for i in range(len(lines_fgin)):
      if not lines_fgin[i]=='':
        word=lines_fgin[i].split()
        if not (word[0]=='#'):
          icount=icount+1
          if icount==1:
            cell = [float(word[0]), float(word[1]), float(word[2]) ]
            cell_duplicate =np.array(cell)
          if icount==2:
             vbox_a= np.array([ float(word[0]),  float(word[1]),  float(word[2]) ])
          if icount==3:
             vbox_b= np.array([ float(word[0]),  float(word[1]),  float(word[2]) ]) 
          if icount==4:
             vbox_c= np.array([ float(word[0]),  float(word[1]),  float(word[2]) ]) 
             #boxini = np.mat (( vbox_a, vbox_b, vbox_c))
             #debug print vbox_a, np.shape(vbox_a)
             boxini = np.vstack ([ vbox_a, vbox_b, vbox_c])

             #                                              (A)
             #Please pay atttention of the fact that boxini=(B), has the same form as the input  file
             #                                              (C)
             #internal box should be box=(A,B,C)=boxini.T
             #debug print boxini.T[:,0]
          if icount==5:
             nat_per_box=int(word[0])
             cell=np.matrix(np.zeros((nat_per_box,3)))
             itype=np.array(np.zeros(nat_per_box))
          if (icount > 5):
             cell[icount-6,0]=float(word[0])
             cell[icount-6,1]=float(word[1])
             cell[icount-6,2]=float(word[2])
             itype[icount-6]=int(word[3])


     #rcell=np.zeros((nat_per_box,3))
     #for i in range(nat_per_box):
     #   for j in range(3):
     #    for k in range(3):
     #     rcell[i,j]=rcell[i,j]+cell[i,k]*boxini[k,j] 
     
     box=boxini.T
     rcell = np.matmul(box, cell.T).T
 
     return nat_per_box, boxini, itype, cell, rcell


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
             print("WARNING: your reper is not right handed.")
             print("WARNING: This is a critical issue. The LAMMPS results are wrong !!!!!") 


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

def set_cell(self, atoms, change=False):
        lammps_cell, self.coord_transform = convert_cell(atoms.get_cell())
        xhi = lammps_cell[0, 0]
        yhi = lammps_cell[1, 1]
        zhi = lammps_cell[2, 2]
        xy = lammps_cell[0, 1]
        xz = lammps_cell[0, 2]
        yz = lammps_cell[1, 2]

        if change:
            cell_cmd = 'change_box all     x final 0 {} y final 0 {} z final 0 {}      xy final {} xz final {} yz final {}'\
                .format(xhi, yhi, zhi, xy, xz, yz)
        else:
            # just in case we'll want to run with a funny shape box, and here command will only happen once, and before any calculation
            if self.parameters.create_box:
                self.lmp.command('box tilt large')
            cell_cmd = 'region cell prism    0 {} 0 {} 0 {}     {} {} {}     units box'\
                .format(xhi, yhi, zhi, xy, xz, yz)

        self.lmp.command(cell_cmd)


def set_lammps_pos(self, atoms):
        pos = atoms.get_positions() / unit_convert("distance", self.units)

        # If necessary, transform the positions to new coordinate system
        if self.coord_transform is not None:
            pos = np.dot(self.coord_transform , np.matrix.transpose(pos))
            pos = np.matrix.transpose(pos)

        # Convert ase position matrix to lammps-style position array
        lmp_positions = list(pos.ravel())

        # Convert that lammps-style array into a C object
        lmp_c_positions =\
            (ctypes.c_double * len(lmp_positions))(*lmp_positions)
#        self.lmp.put_coosrds(lmp_c_positions)
        self.lmp.scatter_atoms('x', 1, 3, lmp_c_positions)


def print_structures_ndm(nat_per_box,nrepeat,boxfin,rcell,cell):

  fout = open("cube_new.gin","w")
  fout.write("%i %i %i\n"%(nrepeat,nrepeat,nrepeat))
  fout.write(" %22.16f  %22.16f  %22.16f\n"%(boxfin[0,0], boxfin[0,1],boxfin[0,2]))
  fout.write(" %22.16f  %22.16f  %22.16f\n"%(boxfin[1,0], boxfin[1,1],boxfin[1,2]))
  fout.write(" %22.16f  %22.16f  %22.16f\n"%(boxfin[2,0], boxfin[2,1],boxfin[2,2]))
  fout.write("%i\n"%(nat_per_box))


  for i in range(nat_per_box):
    fout.write("%20.12f %20.12f %20.12f  1\n"%(cell[i,0],cell[i,1],cell[i,2]))
  fout.close()

  #writing name.in
  fout = open("name_new.in","w")
  fout.write("cube")
  fout.close()

def print_structures_lammps(nat_per_box,itype,lammps_cell,rcell,name_of_lammps_file=None):
  xlo=0.0
  ylo=0.0
  zlo=0.0
  xhi = lammps_cell[0, 0]
  yhi = lammps_cell[1, 1]
  zhi = lammps_cell[2, 2]
  xy = lammps_cell[0, 1]
  xz = lammps_cell[0, 2]
  yz = lammps_cell[1, 2]

  if name_of_lammps_file is None:
   fout = open("cube.lmp","w")
  else:
   fout = open(name_of_lammps_file,"w")

  fout.write("# add your comment ... \n")
  fout.write(" \n")
  fout.write("%i atoms\n"%nat_per_box)
  ntypes=len(np.unique(itype))
  fout.write("%i atom types\n"%ntypes)
  fout.write(" \n")
  fout.write(" %22.16f  %22.16f   xlo xhi\n"%(xlo, xhi))
  fout.write(" %22.16f  %22.16f   ylo yhi\n"%(ylo, yhi))
  fout.write(" %22.16f  %22.16f   zlo zhi\n"%(zlo, zhi))
  fout.write(" %22.16f  %22.16f  %22.16f   xy xz yz\n"%(xy, xz,  yz))
  fout.write(" \n")
  fout.write("Atoms\n")
  fout.write(" \n")

  for i in range(nat_per_box):
    fout.write("%8i  %4i  %20.12f %20.12f %20.12f\n"%(i+1, itype[i], rcell[i,0],rcell[i,1],rcell[i,2]))
  fout.close()


def usage():
  print('''  Usage:  python gin_to_lammps.py name_of_the_gin_file   name_of_lammps_file ''') 

if  ( len(sys.argv)==1 or len(sys.argv) > 3) :
    usage()
    exit(0)
elif len(sys.argv)==2:
   name_of_gin_file=str(sys.argv[1])
   if not os.path.exists (name_of_gin_file):
     print("gin file doesn't exists. The program will stop")
     exit(0)
   name_of_lammps_file=input("Give the name of the lammps file: ")
elif len(sys.argv) ==3:
   name_of_gin_file=str(sys.argv[1])
   if not os.path.exists (name_of_gin_file):
     print("gin file doesn't exists. The program will stop")
     exit(0)
   name_of_lammps_file=str(sys.argv[2])


nat_per_box, boxini, itype, cell, rcell=read_gin('.',name_of_gin_file)



#boxini is in the form that we see the matrix in input:
#                      a1x, a1y, a1z = a1 = A
#     boxini    =      a2x, a2y, a2z = a2 = B
#                      a3x, a3y, a3z = a3 = C
#
#for the calculus should be used box=boxini.T
#                      a1x, a2x, a3x = 
#     box       =      a1y, a2y, a3y = a1 a2 a3 = A B C
#                      a1z, a2z, a3z = 
#
#box are in the format  box[0-2,:]=(A,B,C) first index cartesian projections, 
#                                          second index axis

#x_cart (3,nat) =  box    (3 ,3) x x_cryst(3,nat)
#x_cryst(3,nat) =  box^-1 (3 ,3) x x_cart (3,nat)




tri_mat, beta_matrix = convert_cell_from_gin_to_lammps(boxini)
# beta_matrix = < e_i', e_j> e_j and e_i' define the reper of old and new axis, repectively.   
# tri_max is in box format

 
transd=np.zeros((3,3))

#for i in range(3):
#  for j in range(3):
#    transd[i,j]=np.dot(tri_mat.T[:,i],boxini.T[:,j])  
#
#these two are equivalent ...
#
#trans=np.matmul(tri_mat,boxini.T)

rcell_new=np.matmul(beta_matrix,rcell.T).T
tri_mat_inv=np.linalg.inv(tri_mat)
cell_new = np.matmul(tri_mat_inv,rcell_new.T).T


#this is just for testing .... 
#box_inv=np.linalg.inv(boxini.T)
#cell_new = np.dot(box_inv,rcell.T).T
#print_structures_ndm(nat_per_box,1,boxini,rcell,cell_new)

print_structures_ndm(nat_per_box,1,tri_mat.T,rcell_new,cell_new)

print_structures_lammps(nat_per_box, itype,  tri_mat, rcell_new,name_of_lammps_file)

