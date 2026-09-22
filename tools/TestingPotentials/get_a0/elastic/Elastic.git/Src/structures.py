import globalv
import numpy as np
import os
from ndm    import read_gin, print_structures_ndm
from vasp   import print_structures_vasp
from pwscf  import print_structures_pwscf
from phondy import print_structures_phondy
from lammps import *
import numpy.matlib

def gen_structure(box,strain_tensor, ginFile=None,outFile=None,generic=None):
# gen_structure_prim is not tested in this new version with lammps

  structure=globalv.structure
  mode=globalv.mode
  if (mode == 'vasp') or (mode == 'bvasp'):

    nat_per_box, boxfin, itype, rcell, cell = gen_structure_cube(box,strain_tensor,structure,ginFile)
    #tested for bcc, fcc, c15
    #nat_per_box, boxfin, itype, rcell, cell = gen_structure_prim(box,strain_tensor,structure)
    nrepeat=4
    print_structures_vasp(nat_per_box,nrepeat,boxfin,itype, rcell,cell)

  if (mode == 'pwscf') or (mode == 'bpwscf'):
    if structure != 'c15':
      nat_per_box, boxfin, itype, rcell, cell = gen_structure_cube(box,strain_tensor,structure,ginFile)
    #tested for bcc, fcc, c15
    if structure=='c15':
      nat_per_box, boxfin, itype, rcell, cell = gen_structure_prim(box,strain_tensor,structure,ginFile)
    nrepeat=4
    print_structures_pwscf(nat_per_box,nrepeat,boxfin,itype, itype, rcell,cell)


  if (mode == 'ndm') :
    if  ( structure != 'c15'):
       nat_per_box, boxfin, itype, rcell, cell = gen_structure_cube(box,strain_tensor,structure,ginFile)
    if  ( structure == 'c15'):
       nat_per_box, boxfin, itype, rcell, cell = gen_structure_prim(box,strain_tensor,structure)

    #tested for bcc, fcc, c15
    #nat_per_box, boxfin, rcell, cell = gen_structure_prim(box,strain_tensor,structure)
    nrepeat=1
    if ( structure != 'gin' ):
      nrepeat=4
    if structure == 'c15':
        nrepeat=4

    print_structures_ndm(nat_per_box,nrepeat,boxfin,itype,rcell,cell
)
  if (mode == 'lammps') or (mode == 'blammps'):

    if  ( structure != 'c15'):
       nat_per_box, boxfin, itype, rcell, cell = gen_structure_cube(box,strain_tensor,structure,ginFile)
    if  ( structure == 'c15'):
       nat_per_box, boxfin, itype, rcell, cell = gen_structure_prim(box,strain_tensor,structure)


    #boxini and boxfin are in the form that we see the matrix in input:
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

    tri_mat, beta_matrix = convert_cell_from_gin_to_lammps(boxfin)
    # beta_matrix = < e_i', e_j> e_j and e_i' define the reper of old and new axis, repectively.
    # tri_max is in box format
    if beta_matrix is None:
       beta_matrix=np.matlib.eye(3,dtype=float)

    globalv.beta_matrix=beta_matrix

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

    print_structures_lammps(nat_per_box, itype,  tri_mat, rcell_new,outFile)



    #tested for bcc, fcc, c15
    #nat_per_box, boxfin, rcell, cell = gen_structure_prim(box,strain_tensor,structure)
    nrepeat=1
    if ( structure != 'gin' ):
       nrepeat=4
    if structure == 'c15':
       nrepeat=4


  if (mode == 'phondy') or (mode == 'bphondy'):
    if  ( structure != 'c15'):
       nat_per_box, boxfin, itype, rcell, cell = gen_structure_cube(box,strain_tensor,structure,ginFile)
    if  ( structure == 'c15'):
       nat_per_box, boxfin, itype, rcell, cell = gen_structure_prim(box,strain_tensor,structure)

    nrepeat=6
    print_structures_phondy(nat_per_box,nrepeat,boxfin,itype, rcell,cell)

    globalv.atoms_per_unit=nat_per_box
  if generic is None:
    return;
  else:
     return nat_per_box,boxfin,itype, rcell,cell;

def gen_structure_cube(box,strain_tensor,lattice_type,ginFile):
#Generate crystal structure in VASP format.
  ax,ay,az = box[0],box[1],box[2]
  e1,e2,e3,e4,e5,e6 = strain_tensor

  if lattice_type=='c15':
   aini_x = [       ax,       0.0,     0.0]
   aini_y = [      0.0,        ay,     0.0]
   aini_z = [      0.0,       0.0,      az]

   boxini = np.asmatrix( (aini_x,aini_y,aini_z))
   nat_per_box=24
   itype=np.full(nat_per_box,1,dtype=int)
 # the cartesian coordinates in ax units. In all others
 # structures fcc, bcc we have inserted the fractional coordinates ...
   cell = np.asmatrix(( (0.000,     0.000,      0.000),\
      (0.250,     0.250,      0.250),\
      (0.500,     0.500,      0.000),\
      (0.500,     0.000,      0.500),\
      (0.000,     0.500,      0.500),\
      (0.750,     0.750,      0.250),\
      (0.625,     0.125,      0.125),\
      (0.750,     0.250,      0.750),\
      (0.875,     0.375,      0.125),\
      (0.875,     0.125,      0.375),\
      (0.625,     0.375,      0.375),\
      (0.125,     0.625,      0.125),\
      (0.250,     0.750,      0.750),\
      (0.375,     0.875,      0.125),\
      (0.375,     0.625,      0.375),\
      (0.125,     0.875,      0.375),\
      (0.625,     0.625,      0.625),\
      (0.375,     0.375,      0.625),\
      (0.375,     0.125,      0.875),\
      (0.125,     0.375,      0.875),\
      (0.875,     0.875,      0.625),\
      (0.875,     0.625,      0.875),\
      (0.625,     0.875,      0.875),\
      (0.125,     0.125,      0.625) ))
  if lattice_type=='a15':
   aini_x = [       ax,       0.0,     0.0]
   aini_y = [      0.0,        ay,     0.0]
   aini_z = [      0.0,       0.0,      az]

   boxini = np.asmatrix( (aini_x,aini_y,aini_z))
   nat_per_box=8
   itype=np.full(nat_per_box,1,dtype=int)
 # the cartesian coordinates in ax units. In all others
 # structures fcc, bcc we have inserted the fractional coordinates ...
   cell = np.asmatrix(( (0.000,     0.000,      0.000),\
      (0.500,     0.500,      0.500),\
      (0.250,     0.500,      0.000),\
      (0.750,     0.500,      0.000),\
      (0.250,     0.000,      0.500),\
      (0.000,     0.750,      0.500),\
      (0.500,     0.000,      0.250),\
      (0.500,     0.000,      0.750) ))

  if  lattice_type == 'bcc':
   aini_x = [       ax,       0.0,     0.0]
   aini_y = [      0.0,        ay,     0.0]
   aini_z = [      0.0,       0.0,      az]

   boxini = np.asmatrix( (aini_x,aini_y,aini_z))
   nat_per_box=2
   itype=np.full(nat_per_box,1,dtype=int)

   cell = np.asmatrix(( ( 0.000,    0.000,   0.000),\
                   ( 0.500,    0.500,   0.500)))
  if  lattice_type == 'fcc':
   aini_x = [       ax,       0.0,     0.0]
   aini_y = [      0.0,        ay,     0.0]
   aini_z = [      0.0,       0.0,      az]

   boxini = np.asmatrix( (aini_x,aini_y,aini_z))
   nat_per_box=4
   itype=np.full(nat_per_box,1,dtype=int)
   cell = np.asmatrix(( ( 0.000,    0.000,   0.000),\
                   ( 0.500,    0.500,   0.000),\
                   ( 0.000,    0.500,   0.500),\
                   ( 0.500,    0.000,   0.500)))


  if  lattice_type == 'gin':
   nat_per_box, boxini, itype, cell, rcell=read_gin(ginFile=ginFile)
   #old root_dir=globalv.root_dir
   #old dirSTR=root_dir+'/Structure'
   #old if ginFile is None:
   #old   sgin= dirSTR+'/'+ 'structure.gin'
   #old else:
   #old   sgin= dirSTR+'/'+ ginFile

   #old if not os.path.exists(sgin):
   #old    print "Gin file %s doesn't exist. Put the correct path in setup_ndm and structure"%str(sgin)
   #old    exit(0)

   #old # should be readed in gin ...
   #old fgin=open(sgin,'r')
   #old lines_fgin=fgin.read().splitlines()
   #old icount=0
   #old for i in range(len(lines_fgin)):
   #old  if not lines_fgin[i]=='':
   #old    word=lines_fgin[i].split()
   #old    if not (word[0]=='#'):
   #old      icount=icount+1
   #old      if icount==1:
   #old        cell = [float(word[0]), float(word[1]), float(word[2]) ]
   #old        cell_duplicate =np.array(cell)
   #old      if icount==2:
   #old         vbox_a= [ float(word[0]),  float(word[1]),  float(word[2]) ]
   #old      if icount==3:
   #old         vbox_b= [ float(word[0]),  float(word[1]),  float(word[2]) ]
   #old      if icount==4:
   #old         vbox_c= [ float(word[0]),  float(word[1]),  float(word[2]) ]
   #old         boxini = np.mat (( vbox_a, vbox_b, vbox_c))
   #old         #if I rescale should be here ....
   #old      if icount==5:
   #old         nat_per_box=int(word[0])
   #old         cell=np.matrix(np.zeros((nat_per_box,3)))
   #old         itype=np.array(np.zeros(nat_per_box))
   #old      if (icount > 5):
   #old         cell[icount-6,0]=float(word[0])
   #old         cell[icount-6,1]=float(word[1])
   #old         cell[icount-6,2]=float(word[2])
   #old         itype[icount-6]=int(word[3])

  e_p = np.asmatrix ([[(1.0+e1),0.5*e6,0.5*e5],\
                  [0.5*e6,(1.0+e2),0.5*e4],\
                  [0.5*e5,0.5*e4,(1.0+e3)]])

  #boxfin = np.matrix(np.zeros((3,3)))
  #for i in range(3):
  #  for j in range(3):
  #    for k in range(3):
  #      boxfin[i,j]=boxfin[i,j]+boxini[i,k]*e_p[k,j]
  boxfin=np.matmul(boxini,e_p)
  #print 'boxfin', boxfin [:,1]
  #print np.shape(boxini), np.shape(boxfin)
  rcell=np.zeros((nat_per_box,3))
  #for i in range(nat_per_box):
  #   for j in range(3):
  #    for k in range(3):
  #     rcell[i,j]=rcell[i,j]+cell[i,k]*boxfin[k,j]
  box=boxfin.T
  rcell = np.matmul(box, cell.T).T



  return nat_per_box, boxfin, itype, rcell, cell




def gen_structure_prim(box,strain_tensor,lattice_type):
#Generate crystal structure in VASP format.

  ax,ay,az = box[0],box[1],box[2]
  e1,e2,e3,e4,e5,e6 = strain_tensor

  if lattice_type=='c15':
   aini_x = [    0.0,  ax/2.0, ax/2.0]
   aini_y = [ ay/2.0,     0.0, ay/2.0]
   aini_z = [ az/2.0,  az/2.0,    0.0]

   boxini = np.asmatrix( (aini_x,aini_y,aini_z))
   nat_per_box=6
   itype=np.full(nat_per_box,1,dtype=int)
 # the cartesian coordinates in ax units. In all others
 # structures fcc, bcc we have inserted the fractional coordinates ...
   inicell = np.asmatrix(( ( 0.125,    0.125,   0.125),\
                     ( 0.875,    0.875,   0.875),\
                     ( 0.500,    0.500,   0.500),\
                     ( 0.500,    0.250,   0.250),\
                     ( 0.250,    0.500,   0.250),\
                     ( 0.250,    0.250,   0.500)) )
   inicell[:,:] = ax*inicell[:,:]
   cell = np.matrix(np.zeros((nat_per_box,3)))

# the inverse of the matrix boxini
   boxiniI =np.matrix(np.zeros((3,3)))
   boxiniI= boxini.I
# the initial coordinates in reduced units N,3 = N,3 x 3,3
   cell = inicell*boxiniI
  if  lattice_type == 'bcc':
   aini_x = [  -ax/2.0,    ax/2.0,  ax/2.0]
   aini_y = [   ay/2.0,   -ay/2.0,  ay/2.0]
   aini_z = [   az/2.0,    az/2.0, -az/2.0]

   boxini = np.asmatrix( (aini_x,aini_y,aini_z))
   nat_per_box=1
   itype=np.full(nat_per_box,1,dtype=int)
   cell = np.asmatrix(( ( 0.000,    0.000,   0.000)) )
  if  lattice_type == 'fcc':
   aini_x = [      0.0,    ax/2.0,  ax/2.0]
   aini_y = [   ay/2.0,       0.0,  ay/2.0]
   aini_z = [   az/2.0,    az/2.0,     0.0]

   boxini = np.asmatrix( (aini_x,aini_y,aini_z))
   nat_per_box=1
   cell = np.asmatrix(( ( 0.000,    0.000,   0.000)) )


  e_p = np.asmatrix ([[(1.0+e1),0.5*e6,0.5*e5],\
                  [0.5*e6,(1.0+e2),0.5*e4],\
                  [0.5*e5,0.5*e4,(1.0+e3)]])

  boxfin = np.matrix(np.zeros((3,3)))
  for i in range(3):
    for j in range(3):
      for k in range(3):
        boxfin[i,j]=boxfin[i,j]+boxini[i,k]*e_p[k,j]
  rcell=np.zeros((nat_per_box,3))
  for i in range(nat_per_box):
     for j in range(3):
      for k in range(3):
       rcell[i,j]=rcell[i,j]+cell[i,k]*boxfin[k,j]
  return nat_per_box, boxfin, itype, rcell, cell

def print_structures_phondy(nat_per_box,nrepeat,boxfin,itype,rcell,cell,outFile=None):
  if np.amax(itype) > 1:
    print('WARNING: print_structure_phondy doesnt handle more than one type')
    exit(0)

  if outFile is None:
    fout = open("geom.data","w")
  else:
    fout = open(outFile,"w")

  fout.write("%i\n"%(nat_per_box))
  fout.write("1.0\n")
  fout.write(" %22.16f  %22.16f  %22.16f\n"%(boxfin[0,0], boxfin[0,1],boxfin[0,2]))
  fout.write(" %22.16f  %22.16f  %22.16f\n"%(boxfin[1,0], boxfin[1,1],boxfin[1,2]))
  fout.write(" %22.16f  %22.16f  %22.16f\n"%(boxfin[2,0], boxfin[2,1],boxfin[2,2]))

  for i in range(nat_per_box):
    fout.write("%15.7f %15.7f %15.7f\n"%(rcell[i,0],rcell[i,1],rcell[i,2]))

  fout.close()

def print_structures_ndm(nat_per_box,nrepeat,boxfin,itype,rcell,cell):

  fout = open("cube.gin","w")
  fout.write("%i %i %i\n"%(nrepeat,nrepeat,nrepeat))
  fout.write(" %22.16f  %22.16f  %22.16f\n"%(boxfin[0,0], boxfin[0,1],boxfin[0,2]))
  fout.write(" %22.16f  %22.16f  %22.16f\n"%(boxfin[1,0], boxfin[1,1],boxfin[1,2]))
  fout.write(" %22.16f  %22.16f  %22.16f\n"%(boxfin[2,0], boxfin[2,1],boxfin[2,2]))
  fout.write("%i\n"%(nat_per_box))


  for i in range(nat_per_box):
    fout.write("%15.7f %15.7f %15.7f  %4i\n"%(cell[i,0],cell[i,1],cell[i,2],itype[i]))
  fout.close()

  #writing name.in
  fout = open("name.in","w")
  fout.write("cube")
  fout.close()

def print_structures_vasp(nat_per_box,nrepeat,boxfin,itype,rcell,cell):
  if np.amax(itype) > 2:
    print('WARNING: print_structure_vasp doesnt handle more than two types')
    exit(0)

  fout = open("POSCAR","w")
  fout.write("W\n")
  fout.write("1.0\n")
  fout.write(" %22.16f  %22.16f  %22.16f\n"%(boxfin[0,0], boxfin[0,1],boxfin[0,2]))
  fout.write(" %22.16f  %22.16f  %22.16f\n"%(boxfin[1,0], boxfin[1,1],boxfin[1,2]))
  fout.write(" %22.16f  %22.16f  %22.16f\n"%(boxfin[2,0], boxfin[2,1],boxfin[2,2]))

  if np.amax(itype) == 2:
   fout.write("%i %i\n"%(np.count_nonzero(itype == 1), np.count_nonzero(itype ==2)))
  if np.amax(itype) == 1:
   fout.write("%i\n"%(nat_per_box))
  fout.write("Direct\n")
  for i in range(nat_per_box):
    fout.write("%15.7f %15.7f %15.7f   \n"%(cell[i,0],cell[i,1],cell[i,2]))
  fout.close()

def print_structures_pwscf(nat_per_box,nrepeat,boxfin,itype,rcell,cell):
  if np.amax(itype) > 1:
    print('WARNING: print_structure_pwscf doent handle more than one type')
    exit(0)
  fout = open("POSCAR.pwscf","w")
  fout.write("CELL_PARAMETERS { angstrom } \n")
  fout.write(" %22.16f  %22.16f  %22.16f\n"%(boxfin[0,0], boxfin[0,1],boxfin[0,2]))
  fout.write(" %22.16f  %22.16f  %22.16f\n"%(boxfin[1,0], boxfin[1,1],boxfin[1,2]))
  fout.write(" %22.16f  %22.16f  %22.16f\n"%(boxfin[2,0], boxfin[2,1],boxfin[2,2]))
  #fout.write("%i\n"%(nat_per_box))
  fout.write("ATOMIC_POSITIONS { crystal }\n")
  for i in range(nat_per_box):
    fout.write("Fe   %15.7f %15.7f %15.7f   \n"%(cell[i,0],cell[i,1],cell[i,2]))
  fout.close()

def print_structures_lammps(nat_per_box,itype,lammps_cell,rcell,NameFile=None):
  xlo=0.0
  ylo=0.0
  zlo=0.0
  xhi = lammps_cell[0, 0]
  yhi = lammps_cell[1, 1]
  zhi = lammps_cell[2, 2]
  xy = lammps_cell[0, 1]
  xz = lammps_cell[0, 2]
  yz = lammps_cell[1, 2]

  if NameFile is None:
    fout = open("cube.lmp","w")
  else:
    fout=open(NameFile,"w")

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
